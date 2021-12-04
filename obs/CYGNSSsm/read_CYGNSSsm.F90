!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_CYGNSSsm
! \label{read_CYGNSSsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5.
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data.
!  15 Jul 2018: Mahdi Navari: Bug in the SMAP reader was fixed
!  31 Aug 2018: Mahdi Navari, Edited to read SPL3SMP.005 & SPL3SMP_E.002
!  1  Apr 2019: Yonghwan Kwon: Upated for reading monthy CDF for the current month
! 11 July 2019: Mahdi Navari, There are several version of SMAP sm data available in each directory
!                  with different Release number and different CRID Version Number. The reader was
!                  modified to read the latest version of data (the reader no longer reads the symbolic
!                  link to the SMAP sm data).
!  8 July 2020: David Mocko: Removed config entry to toggle the QC check.
!                            The QC is now always ON for NASA SMAP SM DA.
!  11 Aug 2020: Yonghwan Kwon: Incorporated Sujay's modifications to support SMAP L2 assimilation
!  28 Jan 2021: Hyunglok Kim: Modified to use CYGNS daily/subdaily/SMAP-combined DA
!
! !INTERFACE:
subroutine read_CYGNSSsm(n, k, OBS_State, OBS_Pert_State)
! !USES:
   use ESMF
   use LIS_mpiMod
   use LIS_coreMod
   use LIS_logMod
   use LIS_timeMgrMod
   use LIS_dataAssimMod
   use LIS_DAobservationsMod
   use map_utils
   use LIS_pluginIndices
   use CYGNSSsm_Mod, only: CYGNSSsm_struc

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: k
   type(ESMF_State)    :: OBS_State
   type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!
!  reads the AMSRE soil moisture observations
!  from NETCDF files and applies the spatial masking for dense
!  vegetation, rain and RFI. The data is then rescaled
!  to the land surface model's climatology using rescaling
!  algorithms.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
   real, parameter        ::  minssdev = 0.05
   real, parameter        ::  maxssdev = 0.11
   real, parameter       :: MAX_SM_VALUE = 0.45, MIN_SM_VALUE = 0.0001
   integer                :: status
   integer                :: grid_index
   character*100          :: smobsdir
   character*100          :: fname
   logical                :: alarmCheck, file_exists
   integer                :: t, c, r, i, j, p, jj
   real,          pointer :: obsl(:)
   type(ESMF_Field)       :: smfield, pertField
   integer                :: gid(LIS_rc%obs_ngrid(k))
   integer                :: assimflag(LIS_rc%obs_ngrid(k))
   real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
   logical                :: data_update
   logical                :: data_upd_flag(LIS_npes)
   logical                :: data_upd_flag_local
   logical                :: data_upd
   real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
   real                   :: smobs_daily(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
   real                   :: sm_current(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))
   real                   :: dt
   real                   :: lon
   real                   :: lhour
   real                   :: gmt
   integer                :: zone
   integer                :: fnd
   real, allocatable      :: ssdev(:)
   integer                :: lis_julss
   real                   :: smvalue
   real                   :: model_delta(LIS_rc%obs_ngrid(k))
   real                   :: obs_delta(LIS_rc%obs_ngrid(k))
   character*4            :: yyyy
   character*8            :: yyyymmdd
   character*2            :: mm, dd, hh
   integer                :: yr, mo, da, hr, mn, ss
   integer                :: cyr, cmo, cda, chr, cmn, css
   integer                :: nyr, nmo, nda, nhr, nmn, nss
   real*8                 :: timenow, time1,time2,time3
   integer                :: doy
   character*200          :: list_files
   character*100          :: temp1
   character*1            :: fproc(4)
   character(len=4)       :: istring
   character*100          :: cyg_subdaily_filename(10)
   character*100          :: cygsmap_filename(10)
   integer                :: mn_ind
   integer                :: ftn, ierr
   character(len=200)     :: cmd
   integer                :: rc
   character(len=3)       :: CRID
   integer, external      :: create_filelist ! C function


   call ESMF_AttributeGet(OBS_State, "Data Directory", &
                          smobsdir, rc=status)
   call LIS_verify(status)
   call ESMF_AttributeGet(OBS_State, "Data Update Status", &
                          data_update, rc=status)
   call LIS_verify(status)

   data_upd = .false.
   obs_unsc = LIS_rc%udef

   alarmCheck = LIS_isAlarmRinging(LIS_rc, "CYGNSS read alarm")

   smobs_daily = LIS_rc%udef

   cyr = LIS_rc%yr
   cmo = LIS_rc%mo
   cda = LIS_rc%da
   chr = LIS_rc%hr
   cmn = LIS_rc%mn
   css = LIS_rc%ss

   call LIS_tick(time1,doy,gmt,cyr,cmo,cda,chr,cmn,css,0.0)
   nyr = LIS_rc%yr
   nmo = LIS_rc%mo
   nda = LIS_rc%da
   nhr = LIS_rc%hr
   nmn = LIS_rc%mn
   nss = LIS_rc%ss

   call LIS_tick(time2,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,3600.0)
   nyr = LIS_rc%yr
   nmo = LIS_rc%mo
   nda = LIS_rc%da
   nhr = LIS_rc%hr
   nmn = LIS_rc%mn
   nss = LIS_rc%ss

   call LIS_tick(time3,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,LIS_rc%ts)

   if (alarmCheck .or. CYGNSSsm_struc(n)%startMode) then
      CYGNSSsm_struc(n)%startMode = .false.
      if  (CYGNSSsm_struc(n)%data_designation.eq."CYG_SubDaily") then 

         CYGNSSsm_struc(n)%smobs = LIS_rc%udef
         CYGNSSsm_struc(n)%smtime = -1.0
 
         write(temp1,fmt='(i4.4)') LIS_localPet
         read(temp1,fmt='(4a1)') fproc
         write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
         write(yyyy,'(i4.4)') LIS_rc%yr
         write(mm,'(i2.2)') LIS_rc%mo
         write(dd,'(i2.2)') LIS_rc%da
         write(hh,'(i2.2)') LIS_rc%hr

         if(LIS_masterproc) then
            list_files = "ls "//trim((smobsdir))//&
                 "/"//trim(yyyy)//"."//trim(mm)//"."//dd//&
                 "/CYG_SubDaily_"//trim(yyyymmdd)//"T"//trim(hh)&
                 //"*.h5 > CYGNSS_SubDaily_filelist_sm.dat"

            call system(trim(list_files))
            do i=0,LIS_npes-1
               write(istring,'(I4.4)') i
               cmd = 'cp CYGNSS_SubDaily_filelist_sm.dat CYGNSS_SubDaily_filelist.sm.'//istring//".dat"
               call system(trim(cmd))
            end do ! i
         end if
#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm,ierr)
#endif

         i = 1
         ftn = LIS_getNextUnitNumber()
         open(ftn,file="./CYGNSS_SubDaily_filelist.sm."//&
              fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat",&
              status='old',iostat=ierr)

         do while(ierr.eq.0)
            read(ftn,'(a)',iostat=ierr) fname
            if(ierr.ne.0) then
               exit
            endif
            !mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))

            mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11
            read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
            ss=0
            call LIS_tick(timenow,doy,gmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                 LIS_rc%hr, mn, ss, 0.0)

            cyg_subdaily_filename(i) = fname

            write(LIS_logunit,*) '[INFO] reading ',trim(cyg_subdaily_filename(i))
          
            call read_CYGNSS_SubDailysm_data(n,k,cyg_subdaily_filename(i),&
                 CYGNSSsm_struc(n)%smobs,timenow)

            i = i+1
         enddo
         call LIS_releaseUnitNumber(ftn)


      elseif (CYGNSSsm_struc(n)%data_designation.eq."CYG_SMAP") then 

         CYGNSSsm_struc(n)%smobs = LIS_rc%udef
         CYGNSSsm_struc(n)%smtime = -1.0
 
         write(temp1,fmt='(i4.4)') LIS_localPet
         read(temp1,fmt='(4a1)') fproc
         write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
         write(yyyy,'(i4.4)') LIS_rc%yr
         write(mm,'(i2.2)') LIS_rc%mo
         write(dd,'(i2.2)') LIS_rc%da
         write(hh,'(i2.2)') LIS_rc%hr

         if(LIS_masterproc) then
            list_files = "ls "//trim((smobsdir))//&
                 "/"//trim(yyyy)//"."//trim(mm)//"."//dd//&
                 "/*_"//trim(yyyymmdd)//"T"//trim(hh)&
                 //"*.h5 > CYGSMAP_filelist_sm.dat"

            call system(trim(list_files))
            do i=0,LIS_npes-1
               write(istring,'(I4.4)') i
               cmd = 'cp CYGSMAP_filelist_sm.dat CYGSMAP_filelist.sm.'//istring//".dat"
               call system(trim(cmd))
            end do ! i
         end if
#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm,ierr)
#endif

         i = 1
         ftn = LIS_getNextUnitNumber()
         open(ftn,file="./CYGSMAP_filelist.sm."//&
              fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat",&
              status='old',iostat=ierr)

         do while(ierr.eq.0)
            read(ftn,'(a)',iostat=ierr) fname
            if(ierr.ne.0) then
               exit
            endif

            mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11
            read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
            ss=0
            call LIS_tick(timenow,doy,gmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                 LIS_rc%hr, mn, ss, 0.0)

            cygsmap_filename(i) = fname

            write(LIS_logunit,*) '[INFO] reading ',trim(cygsmap_filename(i))
         
            if (index(cygsmap_filename(i), 'CYG').gt.0) then
            call read_CYGNSS_SubDailysm_data(n,k,cygsmap_filename(i),&
                 CYGNSSsm_struc(n)%smobs,timenow)
            
            elseif (index(cygsmap_filename(i), 'SMAP_L2').gt.0) then
            call read_SMAPL2sm_data_in_CYG(n,k,cygsmap_filename(i),&
                 CYGNSSsm_struc(n)%smobs,timenow)
        
            endif

            i = i+1
         enddo
         call LIS_releaseUnitNumber(ftn)


      elseif (CYGNSSsm_struc(n)%data_designation .eq. "CYG_Daily") then
!-------------------------------------------------------------------------
! MN: create filename for 36km product  (CYGNSS Daily)
!-------------------------------------------------------------------------

         write (yyyy, '(i4.4)') LIS_rc%yr
         write (mm, '(i2.2)') LIS_rc%mo
         write (dd, '(i2.2)') LIS_rc%da

         if (LIS_masterproc) then
            list_files = trim(smobsdir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                         trim(dd)//'/CYG_Daily_' &
                         //trim(yyyy)//trim(mm)//trim(dd)&
                         //'*.h5'
            write(LIS_logunit,*) &
                  '[INFO] Searching for ',trim(list_files)
            rc = create_filelist(trim(list_files)//char(0), &
                 "CYGNSS_Daily_filelist.dat"//char(0))
            if (rc .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] Problem encountered when searching for SMAP files'
               write(LIS_logunit,*) &
                    'Was searching for ',trim(list_files)
               write(LIS_logunit,*) &
                    'LIS will continue...'
            endif
         end if
#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm, ierr)
#endif

         ftn = LIS_getNextUnitNumber()
         open (ftn, file="./CYGNSS_Daily_filelist.dat", &
               action='read', status='old', iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur.
! This assumes that the 'ls command' will list the files in that order.

         do while (ierr .eq. 0)
            read (ftn, '(a)', iostat=ierr) fname
            if (ierr .ne. 0) then
               exit
            endif

            write (LIS_logunit, *) '[INFO] Reading daily CYNGSS ', trim(fname)
            call read_CYGNSS_Dailysm_data(n, k, fname, smobs_daily)

         enddo

         CYGNSSsm_struc(n)%smobs = LIS_rc%udef
         CYGNSSsm_struc(n)%smtime = -1
         call LIS_releaseUnitNumber(ftn)
!-------------------------------------------------------------------------
!   CYGNSS daily data assumed to overpass at 12pm localtime 
!-------------------------------------------------------------------------
         do r = 1, LIS_rc%obs_lnr(k)
            do c = 1, LIS_rc%obs_lnc(k)
               grid_index = LIS_obs_domain(n, k)%gindex(c, r)
               if (grid_index .ne. -1) then

                  if (smobs_daily(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0) then
                     CYGNSSsm_struc(n)%smobs(c, r) = &
                        smobs_daily(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 12.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     CYGNSSsm_struc(n)%smtime(c, r) = gmt

                  endif
               endif
            enddo
         enddo

#if 0
!-------------------------------------------------------------------------
!  From the SMAP documentation:
!  The current approach for the SPL3SMP product is to use the nearest
!  6:00 a.m. LST criterion to perform Level-3 compositing for the
!  descending passes. According to this criterion, for a given grid cell,
!  an L2 data point acquired closest to 6:00 a.m. local solar time will
!  make its way to the final Level-3 file; other late-coming L2 data
!  points falling into the same grid cell will be ignored. For a given
!  file whose time stamp (yyyy-mm-ddThh:mm:ss) is expressed in UTC, only
!  the hh:mm:ss part is converted into local solar time.
!  (O'Neill et al. 2012)
!-------------------------------------------------------------------------
         do r = 1, LIS_rc%obs_lnr(k)
            do c = 1, LIS_rc%obs_lnc(k)
               grid_index = LIS_obs_domain(n, k)%gindex(c, r)
               if (grid_index .ne. -1) then

                  if (smobs(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0) then
                     CYGNSSsm_struc(n)%smobs(c, r) = &
                        smobs(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 6.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     CYGNSSsm_struc(n)%smtime(c, r) = gmt
                  endif
               endif
            enddo
         enddo
#endif

      endif ! sensor
   endif ! alram

   call ESMF_StateGet(OBS_State, "Observation01", smfield, &
                      rc=status)
   call LIS_verify(status, 'Error: StateGet Observation01')

   call ESMF_FieldGet(smfield, localDE=0, farrayPtr=obsl, rc=status)
   call LIS_verify(status, 'Error: FieldGet')

   fnd = 0
   sm_current = LIS_rc%udef

   ! dt is not defined as absolute value of the time difference to avoid
   ! double counting of the data in assimilation.

   if ( (CYGNSSsm_struc(n)%data_designation.eq."CYG_SubDaily") .or. &
        (CYGNSSsm_struc(n)%data_designation.eq."CYG_SMAP") ) then

      do r=1,LIS_rc%obs_lnr(k)
         do c=1,LIS_rc%obs_lnc(k)
            if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
               grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
               dt = (CYGNSSsm_struc(n)%smtime(c,r)-time1)
               if(dt.ge.0.and.dt.lt.(time3-time1)) then
                  sm_current(c,r) = &
                       CYGNSSsm_struc(n)%smobs(c,r)
                  if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                     obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                          sm_current(c,r)
                  endif
                  if(sm_current(c,r).ne.LIS_rc%udef) then
                     fnd = 1
                  endif
               endif
            endif
         enddo
      enddo

   else

      do r = 1, LIS_rc%obs_lnr(k)
         do c = 1, LIS_rc%obs_lnc(k)
            if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
               grid_index = c + (r - 1)*LIS_rc%obs_lnc(k)

               dt = (LIS_rc%gmt - CYGNSSsm_struc(n)%smtime(c, r))*3600.0
               if (dt .ge. 0 .and. dt .lt. LIS_rc%ts) then
                  sm_current(c, r) = &
                        CYGNSSsm_struc(n)%smobs(c, r)
                  if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
                     obs_unsc(LIS_obs_domain(n, k)%gindex(c, r)) = &
                        sm_current(c, r)
                  endif
                  if (sm_current(c, r) .ne. LIS_rc%udef) then
                     fnd = 1
                  endif
               endif
            endif
         enddo
      enddo
   endif

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------

   ! Read monthly CDF (only for the current month)
   if (CYGNSSsm_struc(n)%ntimes .gt. 1 .and. CYGNSSsm_struc(n)%cdf_read_opt .eq. 1) then
      if (.not. CYGNSSsm_struc(n)%cdf_read_mon .or. LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
          LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
         call LIS_readMeanSigmaData(n, k, &
                                    CYGNSSsm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    CYGNSSsm_struc(n)%modelcdffile, &
                                    "SoilMoist", &
                                    CYGNSSsm_struc(n)%model_mu, &
                                    CYGNSSsm_struc(n)%model_sigma, &
                                    LIS_rc%mo)

         call LIS_readMeanSigmaData(n, k, &
                                    CYGNSSsm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    CYGNSSsm_struc(n)%obscdffile, &
                                    "SoilMoist", &
                                    CYGNSSsm_struc(n)%obs_mu, &
                                    CYGNSSsm_struc(n)%obs_sigma, &
                                    LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              CYGNSSsm_struc(n)%nbins, &
                              CYGNSSsm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              CYGNSSsm_struc(n)%modelcdffile, &
                              "SoilMoist", &
                              CYGNSSsm_struc(n)%model_xrange, &
                              CYGNSSsm_struc(n)%model_cdf, &
                              LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              CYGNSSsm_struc(n)%nbins, &
                              CYGNSSsm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              CYGNSSsm_struc(n)%obscdffile, &
                              "SoilMoist", &
                              CYGNSSsm_struc(n)%obs_xrange, &
                              CYGNSSsm_struc(n)%obs_cdf, &
                              LIS_rc%mo)

         CYGNSSsm_struc(n)%cdf_read_mon = .true.
      endif
   endif

   if (LIS_rc%dascaloption(k) .eq. "CDF matching" .and. fnd .ne. 0) then
      if (CYGNSSsm_struc(n)%ntimes .gt. 1 .and. CYGNSSsm_struc(n)%cdf_read_opt .eq. 1) then
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            CYGNSSsm_struc(n)%nbins, &
            1, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            CYGNSSsm_struc(n)%model_xrange, &
            CYGNSSsm_struc(n)%obs_xrange, &
            CYGNSSsm_struc(n)%model_cdf, &
            CYGNSSsm_struc(n)%obs_cdf, &
            sm_current)
      else
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            CYGNSSsm_struc(n)%nbins, &
            CYGNSSsm_struc(n)%ntimes, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            CYGNSSsm_struc(n)%model_xrange, &
            CYGNSSsm_struc(n)%obs_xrange, &
            CYGNSSsm_struc(n)%model_cdf, &
            CYGNSSsm_struc(n)%obs_cdf, &
            sm_current)
      endif
   endif

   obsl = LIS_rc%udef
   do r = 1, LIS_rc%obs_lnr(k)
      do c = 1, LIS_rc%obs_lnc(k)
         if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
            obsl(LIS_obs_domain(n, k)%gindex(c, r)) = sm_current(c, r)
         endif
      enddo
   enddo
   !-------------------------------------------------------------------------
   !  Apply LSM based QC and screening of observations
   !-------------------------------------------------------------------------
   call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+" &
                        //trim(LIS_CYGNSSsmobsId)//char(0), n, k, OBS_state)

   call LIS_checkForValidObs(n, k, obsl, fnd, sm_current)

   if (fnd .eq. 0) then
      data_upd_flag_local = .false.
   else
      data_upd_flag_local = .true.
   endif

#if (defined SPMD)
   call MPI_ALLGATHER(data_upd_flag_local, 1, &
                      MPI_LOGICAL, data_upd_flag(:), &
                      1, MPI_LOGICAL, LIS_mpi_comm, status)
   data_upd = any(data_upd_flag)
#else
   data_upd = data_upd_flag_local
#endif

   if (data_upd) then

      do t = 1, LIS_rc%obs_ngrid(k)
         gid(t) = t
         if (obsl(t) .ne. -9999.0) then
            assimflag(t) = 1
         else
            assimflag(t) = 0
         endif
      enddo

      call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                             .true., rc=status)
      call LIS_verify(status)

      if (LIS_rc%obs_ngrid(k) .gt. 0) then
         call ESMF_AttributeSet(smField, "Grid Number", &
                                gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(smField, "Assimilation Flag", &
                                assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(smfield, "Unscaled Obs", &
                                obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

      endif

      if (LIS_rc%dascaloption(k) .eq. "CDF matching") then
         if (CYGNSSsm_struc(n)%useSsdevScal .eq. 1) then
            call ESMF_StateGet(OBS_Pert_State, "Observation01", pertfield, &
                               rc=status)
            call LIS_verify(status, 'Error: StateGet Observation01')

            allocate (ssdev(LIS_rc%obs_ngrid(k)))
            ssdev = CYGNSSsm_struc(n)%ssdev_inp

            if (CYGNSSsm_struc(n)%ntimes .eq. 1) then
               jj = 1
            else
               jj = LIS_rc%mo
            endif
            do t = 1, LIS_rc%obs_ngrid(k)
               if (CYGNSSsm_struc(n)%obs_sigma(t, jj) .gt. 0) then
                  ssdev(t) = ssdev(t)*CYGNSSsm_struc(n)%model_sigma(t, jj)/ &
                             CYGNSSsm_struc(n)%obs_sigma(t, jj)
                  if (ssdev(t) .lt. minssdev) then
                     ssdev(t) = minssdev
                  endif
               endif
            enddo

            if (LIS_rc%obs_ngrid(k) .gt. 0) then
               call ESMF_AttributeSet(pertField, "Standard Deviation", &
                                      ssdev, itemCount=LIS_rc%obs_ngrid(k), rc=status)
               call LIS_verify(status)
            endif
            deallocate (ssdev)
         endif
      endif
   else
      call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                             .false., rc=status)
      call LIS_verify(status)
   endif

end subroutine read_CYGNSSsm

!BOP
! 
! !ROUTINE: read_CYGNSS_SubDailysm_data
! \label{read_CYGNSS_SubDailysm_data}
!
! !INTERFACE:
subroutine read_CYGNSS_SubDailysm_data(n, k,fname, smobs_inp, time)
! 
! !USES:   

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use CYGNSSsm_Mod, only : CYGNSSsm_struc

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  character (len=*)        :: fname
  real                     :: smobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id
  real,             allocatable  :: sm_field(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  real                           :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')

  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status)
  call LIS_verify(status, 'Error opening CYGNSS SubDaily file ')

  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LIS_verify(status, 'Error opening SM group in CYGNSS SubDaily file')

  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LIS_verify(status, 'Error opening SM field in CYGNSS SubDaily file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LIS_verify(status, 'Error opening row index field in CYGNSS SubDaily file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LIS_verify(status, 'Error opening column index field in CYGNSS SubDaily file')

  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LIS_verify(status, 'Error in h5dget_space_f: rea CYGNSS SubDaily Obs')

! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status)
  if(status.eq.-1) then
     call LIS_verify(status, 'Error in h5sget_simple_extent_dims_f: readCYGNSS SubDailyObs')
  endif

  allocate(sm_field(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LIS_verify(status, 'Error extracting row index from CYGNSS SubDaily file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LIS_verify(status, 'Error extracting col index from CYGNSS SubDaily file')

  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LIS_verify(status, 'Error extracting SM field from CYGNSS SubDaily file')

  call h5dclose_f(row_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5gclose_f(sm_gr_id,status)
  call LIS_verify(status,'Error in H5GCLOSE call')

  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')

  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  sm_data = LIS_rc%udef
  sm_data_b = .false.

!-------------------------------------------------------------------original
!!grid the data in EASE projection
!  do t=1,maxdims(1)
!     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then
!        sm_data(ease_col(t) + &
!             (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = sm_field(t)
!        if(sm_field(t).ne.-9999.0) then
!           sm_data_b(ease_col(t) + &
!                (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = .true.
!        endif
!     endif
!  enddo

!--------------------------------------------------------------------YK
  do t=1,maxdims(1)
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then
        sm_data(ease_col(t) + &
             (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = sm_field(t)
           if(sm_data(ease_col(t) + &
                   (ease_row(t)-1)*CYGNSSsm_struc(n)%nc).ne.-9999.0) then
                 sm_data_b(ease_col(t) + &
                    (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = .true.
           endif
     endif
  enddo
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:), sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       CYGNSSsm_struc(n)%rlat, CYGNSSsm_struc(n)%rlon,&
       CYGNSSsm_struc(n)%n11, LIS_rc%udef, ios)


  deallocate(sm_field)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))

           CYGNSSsm_struc(n)%smtime(c,r) = &
                time
        endif
     enddo
  enddo
#endif

end subroutine read_CYGNSS_SubDailysm_data

!BOP
! 
! !ROUTINE: read_SMAPL2sm_data_in_CYG
! \label{read_SMAPL2sm_data_in_CYG}
!
! !INTERFACE:
subroutine read_SMAPL2sm_data_in_CYG(n, k,fname, smobs_inp, time)
! 
! !USES:   

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use CYGNSSsm_Mod, only : CYGNSSsm_struc

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  character (len=*)        :: fname
  real                     :: smobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"
  character*100,    parameter    :: sm_qa_name = "retrieval_qual_flag"
!YK
  character*100,    parameter    :: vwc_field_name = "vegetation_water_content"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id, sm_qa_id
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  integer(hid_t)                 :: vwc_field_id ! YK
  real,             allocatable  :: sm_field(:)
  real,             allocatable  :: vwc_field(:)! YK
  integer,          allocatable  :: sm_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  real                           :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')

  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status)
  call LIS_verify(status, 'Error opening SMAP L2 file ')

  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LIS_verify(status, 'Error opening SM group in SMAP L2 file')

  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LIS_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LIS_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LIS_verify(status, 'Error opening column index field in SMAP L2 file')

!YK
  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
  call LIS_verify(status, 'Error opening QA field in SMAP L2 file')

!YK 
  call h5dopen_f(sm_gr_id, vwc_field_name,vwc_field_id, status)
  call LIS_verify(status, 'Error opening Veg water content field in SMAP L2 file')

  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LIS_verify(status, 'Error in h5dget_space_f: reaSMAP L2Obs')

! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status)
  if(status.eq.-1) then
     call LIS_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif

  allocate(sm_field(maxdims(1)))
  allocate(sm_qa(maxdims(1)))    !YK
  allocate(vwc_field(maxdims(1)))    !YK
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LIS_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LIS_verify(status, 'Error extracting col index from SMAP L2 file')

  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP L2 file')

!YK
  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP L2 file')

!YK get the vegetation water content
  call h5dread_f(vwc_field_id, H5T_NATIVE_REAL,vwc_field,dims,status)
  call LIS_verify(status, 'Error extracting Veg water content (AM) field from SMAP L2 file')

!YK
  call h5dclose_f(sm_qa_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

!YK 
  call h5dclose_f(vwc_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5gclose_f(sm_gr_id,status)
  call LIS_verify(status,'Error in H5GCLOSE call')

  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')

  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  sm_data = LIS_rc%udef
  sm_data_b = .false.

!-------------------------------------------------------------------original
!!grid the data in EASE projection
!  do t=1,maxdims(1)
!     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then
!        sm_data(ease_col(t) + &
!             (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = sm_field(t)
!        if(sm_field(t).ne.-9999.0) then
!           sm_data_b(ease_col(t) + &
!                (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = .true.
!        endif
!     endif
!  enddo

!--------------------------------------------------------------------YK
!grid the data in EASE projection
! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell. 
! When retrieval is performed, it contains additional bits to further 
! indicate the exit status and quality of the retrieval. The first bit 
! indicates the recommended quality (0-means retrieval has recommended quality).
  do t=1,maxdims(1)
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then
        sm_data(ease_col(t) + &
             (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = sm_field(t)

        if(vwc_field(t).gt.5) then !YK Aply QC : if VWC > 5 kg/m2
           sm_data(ease_col(t) + &
                (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = LIS_rc%udef
        else
           if(sm_data(ease_col(t) + &
                   (ease_row(t)-1)*CYGNSSsm_struc(n)%nc).ne.-9999.0) then
              if(ibits(sm_qa(t),0,1).eq.0) then
                 sm_data_b(ease_col(t) + &
                    (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = .true.
              else
                 sm_data(ease_col(t) + &
                      (ease_row(t)-1)*CYGNSSsm_struc(n)%nc) = LIS_rc%udef
              endif
           endif
        endif
     endif
  enddo
!-----------------------------------------------------------------------

  !t = 1
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:), sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       CYGNSSsm_struc(n)%rlat, CYGNSSsm_struc(n)%rlon,&
       CYGNSSsm_struc(n)%n11, LIS_rc%udef, ios)


  deallocate(sm_field)
!  deallocate(sm_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))

           CYGNSSsm_struc(n)%smtime(c,r) = &
                time
        endif
     enddo
  enddo
#endif

end subroutine read_SMAPL2sm_data_in_CYG

! MN: the data structure in both 36 km and 9 km products is the same therefore  
!         read_NASASMAP_E_data is similar to read_NASASMAP_data

!BOP
! 
! !ROUTINE: read_CYGNSS_data
! \label{read_CYGNSS_data}
!
! !INTERFACE:
subroutine read_CYGNSS_Dailysm_data(n, k, fname, smobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use CYGNSSsm_Mod, only : CYGNSSsm_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMOS NESDIS binary file and applies the data
!  quality flags to filter the data. !\normalsize
!
!  tb_time_seconds
!  Arithmetic average of the same parameters found in the 
!  fore- and aft-looking groups in the input SPL1CTB granule. 
!  The resulting parameter thus describes the average of UTC 
!  acquisition times of SPL1BTB observations whose boresights 
!  fall within a 36 km EASE-Grid 2.0 cell. The result is then 
!  expressed in J2000 seconds (the number of seconds since 
!  11:58:55.816 on January 1, 2000 UT).
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTNASASMAP AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data_Daily"
  character*100,    parameter    :: sm_field_name = "soil_moisture"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2 ! scaler--> rank = 0 ; 2D array--> rank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id
  real,             allocatable  :: sm_field(:,:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/CYGNSSsm_struc(n)%nc, CYGNSSsm_struc(n)%nr/)
  count_file = (/CYGNSSsm_struc(n)%nc, CYGNSSsm_struc(n)%nr/)
  count_mem  = (/CYGNSSsm_struc(n)%nc, CYGNSSsm_struc(n)%nr/)
  
  allocate(sm_field(CYGNSSsm_struc(n)%nc, CYGNSSsm_struc(n)%nr))
  allocate(dims(2))

  dims(1) = CYGNSSsm_struc(n)%nc
  dims(2) = CYGNSSsm_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening CYGNSS file ')
  
  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LIS_verify(status, 'Error opening SM group in CYGNSS file')
     
  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LIS_verify(status, 'Error opening SM field in CYGNSS file')
     
  call h5dget_space_f(sm_field_id, dataspace, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readCYGNSSObs')
     
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       start=offset_file, count=count_file, hdferr=status)
  call LIS_verify(status, 'Error setting hyperslab dataspace in readCYGNSSObs')
     
  call h5screate_simple_f(memrank,dimsm, memspace, status)
  call LIS_verify(status, 'Error in h5create_simple_f; read_CYGNSSsm')
     
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
       start=offset_mem, count=count_mem, hdferr=status)
  call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_CYGNSSsm')
  
  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status, &
       memspace, dataspace)
  call LIS_verify(status, 'Error extracting SM field from CYGNSSfile')

  call h5dclose_f(sm_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5gclose_f(sm_gr_id,status)
  call LIS_verify(status,'Error in H5GCLOSE call')

  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  sm_data_b = .false. 
  t = 1

  do r=1,CYGNSSsm_struc(n)%nr
     do c=1,CYGNSSsm_struc(n)%nc        
      
        sm_data(t) = sm_field(c,r)

        if(sm_data(t).ne.-9999.0) then 
           sm_data_b(t) = .true.
        endif
        t = t+1
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       CYGNSSsm_struc(n)%nc*CYGNSSsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       CYGNSSsm_struc(n)%rlat, CYGNSSsm_struc(n)%rlon,&
       CYGNSSsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_CYGNSS_Dailysm_data

#if 0 
!BOP
! !ROUTINE: create_CYGNSSsm_filename
! \label{create_CYGNSSsm_filename}
! 
! !INTERFACE: 
subroutine create_CYGNSSsm_filename(ndir, designation, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: designation
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the CYGNSS filename (from NSIDC)
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the CYGNSS soil moisture data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated CYGNSS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/CYG_Daily_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
         '.h5'

end subroutine create_CYGNSSsm_filename

#endif



