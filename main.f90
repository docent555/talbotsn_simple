module ic
   use, intrinsic :: iso_c_binding
   use fourier
   use reflection
   use ifport
   use ifcore

   implicit none

   namelist /param/ nz, period, lz, lx, nx, nth, delta, r0, r1, sigma, gamma, xe, c, &
      in_type, central_mirror, cont, it_todo, ops, lambda, mirrx0, cold, interp! a0_peak!,  g_amp, g_x0, g_x1, intrvl, xcp, alfa, x_out, recount, thout, amp_only

   real(c_double) h, lz, lx, hx, hth, delta, imp_x0, imp_xsp, r0, r1, sigma, &
      xe, gamma, c, c3, kappa, lambda, dva_na_k, norm, mirrx0, x_out!, a0_peak, xcp, alfa
   integer(c_int) ::  nx, nth, nz, ix_out, intrvl, it_todo, it_doiter, in_type, it_made = 0, &
                     it_flag = 0, ie, hours, minutes, seconds!, iimp_x0, iimp_xend
   logical(4) pressed
   character(1) key
   integer(c_int), parameter :: esc = 27
   !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
   logical(c_bool) central_mirror, period, ops, cont, cold, interp

   complex(c_double_complex), allocatable :: a1(:), a0(:), ak1(:), ak0(:), atmp(:), jk1(:), jk0(:), k(:), ex(:), dlt(:), tmp(:), &
                                        aktmp(:), akzl(:), k2(:), akz0(:), a0z0(:), a0z0cut(:), ak_z0(:), ak_zl(:), a_z0(:), a_zl(:)
   real(c_double), allocatable :: th0(:, :), th1(:, :), dthdz(:, :), fk1(:), fk2(:), rhs0(:, :), z(:), x(:), g(:)
   integer(c_int), allocatable :: it(:)

   complex(c_double_complex), parameter :: im1 = (0.0d0, 1.0d0)
   real(c_double) :: start_time, finish_time, calc_time

   !interface
   !    subroutine fn_for_fortran_to_call(ptr) &
   !        bind(c, name='fn_for_fortran_to_call')
   !        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
   !        implicit none
   !        type(c_ptr), intent(in), value :: ptr
   !    end subroutine fn_for_fortran_to_call
   !end interface

contains
   subroutine calc_idx()

      implicit none

      dva_na_k = 4.0d0*pi/lambda

      if (period == .true.) then
         lz = lz*(lx*lx)/lambda
      end if

      c3 = c*c*c
      h = lz/nz
      nz = nz + 1

      hth = 2.0d0*pi/nth
      hx = lx/nx

      x_out = lx/2
      ix_out = int(x_out/hx) + 1
      if (ix_out <= 0) ix_out = 1
      if (ix_out > nx) ix_out = nx
      x_out = hx*(ix_out - 1)

      !kakaya po schetu iteraciya
      if (it_flag == 0) it_made = 0
      it_doiter = it_made + it_todo
   end subroutine

   subroutine calc_theta(th, dthdz)

      implicit none

      real(c_double), intent(inout) :: th(:, :), dthdz(:, :)
      integer, dimension(size(th, 1)) :: i
      integer ix

      i = (/1:size(th, 1)/)

      do ix = 1, 2
         th(:, ix) = hth*(i - 1)
         dthdz(:, ix) = delta
      end do

   end subroutine
end module ic

program elektron2dsi
   use ic; use fourier

   implicit none

   interface
      subroutine init() bind(c, name='init')
      end subroutine init
      subroutine finish() bind(c, name='finish')
      end subroutine finish
      function fk_fn(xe) result(fk_res)
         use, intrinsic :: iso_c_binding, only: c_double
         use ic, only: nx
         real(c_double), dimension(nx) :: fk_res
         real(c_double) xe
      end function fk_fn
   end interface

   !type(c_ptr), intent(in), value :: ptr
   real(c_double) eff(2), &
      sum_eff
   integer(c_int) iz, percent, rec_length, nr, first_it, ith
   character*6 str
   integer i

   !schitaem parametry, vychislyaem a0(x), f(x), k2, z, x i t. d.
   call init()

   !to, chto ne zavisit ot j
   dlt = im1*gamma*k2 + sigma   ! k2 - otricatel'noe
   ex = cdexp(-dlt*h)

   rec_length = c_double_complex*2*nx/4

   call cpu_time(start_time)

   !po novomu - dobavochek 23.03.23
   ak0(:) = fs(a0)

   !do i = 1, nx
   !   fk1(:) = fk_fn(hx*(i - 1))
   !   tmp(i) = mysum(fk1*ak0)*dsqrt(2.0d0/nx)
   !end do
   !open (1, file='test.dat')
   !do i = 1, nx
   !   write (1, '(i,4f17.8)') i, real(a0(i)), real(tmp(i)), imag(a0(i)), imag(tmp(i))
   !end do
   !close (1)
   !stop

   first_it = it_made + 1
   do_iter: do it_made = first_it, it_doiter
      percent = 0
      write (*, '(a11,\,i5,\,a2,i5,a,\)') 'iteration #', it_made, ': ', percent, '%'

      if (cold == .true.) then
         ak1(:) = ak0*cdexp(-dlt*lz)
      else

         !nachal'nye usloviya dlya fazy i eyo rasstrojki pri z=0
         call calc_theta(th0, dthdz)

         !raschet kpd v tochke z=0
         eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !(xcp1)
         eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !(xcp2)

         !nachal'nyj tok (z=0)
         jk0 = fk1*jf(th0(:, 1)) + fk2*jf(th0(:, 2))

         do_z: do iz = 1, nz - 1
            rhs0 = rhs(ak0, th0)
            th1 = th0 + dthdz*h + h/2.0d0*rhs0*h !prediktor theta
            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !prediktor toka

            !predictor
            !call step(ak0, ak1, ex, c3, jk0, jk1, dlt, h)
            call step()
            
            !!prediktor a (interpolyaciya)            
            !ak1(2:nx:2) = ak0(2:nx:2)*ex(2:nx:2) + &
            !              c3*(jk0(2:nx:2) + jk1(2:nx:2)*(-1.0d0 + dlt(2:nx:2)*h) + &
            !                  ex(2:nx:2)*(jk1(2:nx:2) - jk0(2:nx:2)*(1.0d0 + dlt(2:nx:2)*h)))/dlt(2:nx:2)/dlt(2:nx:2)/h

            !!prediktor a (trapecii)
            !atmp = (ak0 + c3*h/2.0d0*jk0)*cdexp(-dlt*h) !chast' a
            !ak1 = atmp + c3*h/2.0d0*jk1 !prediktor a

            !a1 = ifs(ak1) !vozvrashchaemsya v real'nost'

            !korrektor theta - 1
            th1 = th0 + dthdz*h + h/6.0d0*rhs0*h &
                  + h/3.0d0*rhs(ak1, th1)*h

            !theta(:, :, iz + 1) = th1 !zapishem traektoriyu elektrona

            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !korrektor toka

            !to, chto zavisit ot toka
            !jkd = jk1 - jk0

            !!korrektor
            !call step(ak0, ak1, ex, c3, jk0, jk1, dlt, h)
            call step()
            
            !!korrektor a (interpolyaciya)
            !ak1(2:nx:2) = ak0(2:nx:2)*ex(2:nx:2) + &
            !              c3*(jk0(2:nx:2) + jk1(2:nx:2)*(-1.0d0 + dlt(2:nx:2)*h) + &
            !                  ex(2:nx:2)*(jk1(2:nx:2) - jk0(2:nx:2)*(1.0d0 + dlt(2:nx:2)*h)))/dlt(2:nx:2)/dlt(2:nx:2)/h

            !!korrektor a (trapecii)
            !atmp = (ak0 + c3*h/2.0d0*jk0)*cdexp(-dlt*h) !chast' a
            !ak1 = atmp + c3*h/2.0d0*jk1 !korrektor a

            dthdz = dthdz + h/2.0d0*(rhs0 + rhs(ak1, th1)) !korrektor theta - 1 end

            !raschet kpd v tochke z = iz*h
            eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !schitaem kpd (xcp1)
            eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !schitaem kpd (xcp2)

            th0 = th1 !dlya sleduyushchego shaga po z
            jk0 = jk1 !dlya sleduyushchego shaga po z
            ak0 = ak1 !dlya sleduyushchego shaga po z
            !dthdz0 = dthdz1 !dlya sleduyushchego shaga po z

            percent = int(real(iz)/real(nz - 2)*100 + 0.5)
            write (*, '(\,a6,i5,a)') '\b\b\b\b\b\b'c, percent, '%'
         end do do_z

         !raschet summy kpd po x v tochke z=lz
         sum_eff = (eff(1) + eff(2))/2.0d0
      end if

      !vozvrashchaemsya v real'nost'
      a1 = ifs(ak1)
      ak_zl(:) = ak1
      a_zl(:) = a1

      !VARIANT 1      
      call maggot_cmplx(ak1, tmp)
      tmp(:) = tmp(:)*cdexp(-dlt*lz)
      !dlya vyvoda v file
      ak_z0(:) = tmp(:)
      a_z0(:) = ifs(tmp)
      !otrazhenie ot pervogo zerkala
      call maggot_cmplx(tmp, ak0)

      !!VARIANT 2      
      !call maggot_cmplx(ak1, tmp)
      !ak_zl(:) = tmp
      !a_zl(:) = ifs(tmp)
      !tmp(:) = tmp*cdexp(-dlt*lz)
      !!otrazhenie ot pervogo zerkala
      !call maggot_cmplx(tmp, ak0)
      !!dlya vyvoda v file
      !ak_z0(:) = ak0(:)
      !a_z0(:) = ifs(ak0)

      write (*, '(a,\)') char(13)
      write (*, '(a,i6,a,f7.3,a,e17.8,a,e17.8,a)') 'iteration #', it_made, ':     a(', x_out, ') = ', &
         cdabs(a1(ix_out)), '   eff = ', sum_eff

      pressed = peekcharqq()
      if (pressed) then
         key = getcharqq()
         if (ichar(key) .eq. esc) then
            write (*, '(/,a)') 'Quit?'
            key = getcharqq()
            if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
               exit
            end if
         end if
      end if

      !call fn_for_fortran_to_call(ptr)
   end do do_iter

   write (*, '(/,/)')

   call cpu_time(finish_time)
   print *, 'computation time = ', finish_time - start_time, ' seconds'

   call write_result()
   call finish()

   print *, 'calculation finished.'
   pause
   stop
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'

contains
   function jf(th)
      implicit none
      real(c_double), intent(in) :: th(:)
      complex(c_double_complex) jf

      jf = 2.0d0/dble(nth)*sum(cdexp(-im1*th))
   end function jf

   function rhs(ak, th)
      implicit none
      complex(c_double_complex), intent(in) :: ak(:)
      real(c_double), intent(in) :: th(:, :)
      real(c_double), dimension(size(th, 1), size(th, 2)) :: rhs

      rhs(:, 1) = dreal(mysum(fk1*ak)*cdexp(im1*th(:, 1)))*dsqrt(2.0d0/nx)
      rhs(:, 2) = dreal(mysum(fk2*ak)*cdexp(im1*th(:, 2)))*dsqrt(2.0d0/nx)

   end function rhs

   function mysum(a)
      implicit none
      complex(c_double_complex), intent(in) :: a(:)
      complex(c_double_complex) :: mysum
      integer(c_int) i, n

      n = size(a)
      mysum = dcmplx(0)

      do i = n, 1, -1
         mysum = mysum + a(i)
      end do
   end function mysum

   function dmysum(a)
      implicit none
      real(c_double), intent(in) :: a(:)
      real(c_double) :: dmysum
      integer(c_int) i, n

      n = size(a)
      dmysum = 0.0d0

      do i = n, 1, -1
         dmysum = dmysum + a(i)
      end do
   end function dmysum
   !end subroutine calculate_fortran
end program elektron2dsi

subroutine init() bind(c, name='init')
   use, intrinsic :: iso_c_binding
   use ic

   implicit none

   integer i

   interface
      subroutine read_param() bind(c, name='read_param')
      end subroutine read_param
      function a0_fn_stat() result(a0_res)
         use ic
         complex(c_double_complex), dimension(nx) :: a0_res
      end function a0_fn_stat
      function fk_fn(xe) result(fk_res)
         use, intrinsic :: iso_c_binding, only: c_double
         use ic, only: nx
         real(c_double), dimension(nx) :: fk_res
         real(c_double) xe
      end function fk_fn
      function k_fn() result(k)
         use ic, only: nx
         use, intrinsic :: iso_c_binding
         complex(c_double_complex), dimension(2*nx) :: k
      end function k_fn
      function k2_fn() result(k2_res)
         use ic
         complex(c_double_complex), dimension(nx) :: k2_res
      end function k2_fn
      function dn_fn() result(dn_res)
         use ic
         complex(c_double_complex), dimension(nx) :: dn_res
      end function dn_fn
   end interface

   call read_param()
   call calc_idx()
   call allocate_arrays()
   call calc_zxit()
   call sincost_init(nx)
   call fft_init(2*nx)
   call reflection_init(nx, lx, mirrx0)
   call dst_init(nx, lx)

   !nachal'nye usloviya dlya a(z=0)
   if (cont .eq. .true.) then
      open (1, file='a0.bin', form='binary', err=101)
      read (1) a0
      close (1)
   else
      a0 = a0_fn_stat()
      !a0z0 = a0
   end if

   !f_k
   fk1(:) = fk_fn(xe)
   fk2(:) = fk_fn(lx - xe)

   !k**2
   if (ops .eq. .false.) then
      k2 = k2_fn()
   else
      k2 = dn_fn()
   end if

   open (1, file='a0.dat')
   do i = 1, nx
      write (1, '(2e17.8,i10)') (i - 1)*hx, dreal(a0(i))
   end do
   close (1)

   open (1, file='fk_and_k2.dat')
   do i = 1, nx
      write (1, '(2e17.8,i10)') fk1(i), fk2(i), int(k2(i))
   end do
   close (1)

   write (*, *) 'nz = ', nz
   write (*, *) 'h = ', h

   print *, 'lz = ', lz
   print *, 'lx = ', lx
   print *, 'c3 = ', c3

   return
101 stop 'error of file open.'
end subroutine init

subroutine finish() bind(c, name='finish')
   use fourier, only: sincost_destroy, fft_destroy

   call sincost_destroy()
   call fft_destroy()
   call deallocate_arrays()
end subroutine finish

subroutine write_result()
   use ic
   use, intrinsic :: iso_c_binding

   implicit none

   integer i, j

   call cpu_time(start_time)

   it_made = it_made - 1

   open (1, file='aend.dat', err=101)
   do i = 1, nx
      write (1, '(1p7e17.8)', err=103) (i - 1)*hx, cdabs(a_z0(i)), cdabs(a_zl(i)), dreal(a_z0(i)), dimag(a_z0(i)), dreal(a_zl(i)), dimag(a_zl(i))
   end do
   close (1)

   open (1, file='akend.dat', err=101)
   do i = 1, nx
      write (1, '(1p7e17.8)', err=103) (i - 1)*hx, cdabs(ak_z0(i)), cdabs(ak_zl(i)), dreal(ak_z0(i)), dimag(ak_z0(i)), dreal(ak_zl(i)), dimag(ak_zl(i))
   end do
   close (1)

   call cpu_time(finish_time)
   print *, 'writing time = ', finish_time - start_time, ' seconds'

   return
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'
end subroutine write_result

subroutine calc_zxit()
   use ic

   implicit none

   integer i

   do i = 1, nz
      z(i) = (i - 1)*h
   end do
   do i = 1, nx
      x(i) = (i - 1)*hx
   end do
   do i = (it_made + 1), it_doiter
      it(i) = i
   end do

   open (1, file='z.dat', err=101)
   do i = 1, nz
      write (1, *, err=103) z(i)
   end do
   close (1)

   open (1, file='x.dat', err=101)
   do i = 1, nx
      write (1, *, err=103) x(i)
   end do
   close (1)

   open (1, file='k.dat', err=101)
   do i = 1, nx
      write (1, *, err=103) i
   end do
   close (1)

   return
101 stop 'error of file open.'
103 stop 'error of file writing.'
end subroutine calc_zxit

function a0_fn_stat() result(a0_res)
   use ic, only: nx, hx, pi, imp_xsp, imp_x0, in_type, lx, central_mirror!, alfa!, a0_peak!, coeff, iimp_x0,, iimp_xend, xcp
   use fourier

   implicit none

   complex(c_double_complex), dimension(nx) :: a0_res, c
   real(c_double), dimension(nx) :: a0env
   integer i, ix(nx)!, icp
   real(8) xcp, a0_peak, alfa

   if (in_type == 1) then
      !nachal'nye usloviya iz garmonik
      c = 0
      c(2) = 25
      c(4) = 17
      c(6) = 21
      c(8) = 12
      c(10) = 7
      c(12) = 16
      a0_res = ifs(c)
      a0_res = cmplx(dreal(a0_res), 0.0d0)
   elseif (in_type == 2) then
      !spec. nachal'nye usloviya dlya a
      !nachal'nye usloviya dlya a (simmetrichnye impul'sy na krayah)
      alfa = 3.000000000000000E-002
      a0_peak = 0.1
      if (central_mirror == .false.) then
         xcp = lx/4
         !icp = xcp/hx + 1
         !iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)
      else
         xcp = lx/2
         !icp = xcp/hx + 1
         !iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx - 1/)

         a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)

      end if

      !open (1, file='test.dat')
      !do i = 1, nx
      !   write (1, '(4e17.8)') (i - 1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !   !write(1, *) ix(i)
      !end do
      !close (1)
      !stop
   else
      print *, 'error: wrong in_type'
      pause
      stop
   end if

end function a0_fn_stat

function fk_fn(xe) result(fk_res)
   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use ic, only: nx, pi, lx

   implicit none

   real(c_double) :: fk_res(nx), xe
   integer(c_int) n(nx), i

   n = (/0:nx - 1/)

   if (xe > 0.0d0) then
      fk_res(:) = dsin(pi*n*xe/lx)
   else
      fk_res(:) = 0.0d0
      print *, 'fk = 0'
   end if

   !open(1, file = 'test.dat')
   !do i=1,nx
   !    write(1,'(i,e17.8,i)') i-1, fk_res(i), n(i)
   !enddo
   !close(1)
   !stop
end function fk_fn

function k_fn() result(k)
   use ic, only: nx
   use, intrinsic :: iso_c_binding

   implicit none

   complex(c_double_complex), dimension(2*nx) :: k
   complex(c_double_complex) :: im1 = (0.0d0, 1.0d0)
   integer nn

   nn = 2*nx

   k = im1*(/0:nn/2 - 1, -nn/2:-1/)
end function k_fn

function k2_fn() result(k2_res)
   use ic

   implicit none

   complex(c_double_complex), dimension(nx) :: k2_res
   integer i
   real(c_double) w

   !k**2
   do i = 1, nx
      w = pi*(i - 1)/lx
      !k2_res(i) = w * w - bylo
      !k2_res(i) = - w * w ! stalo
      k2_res(i) = -w*w/dva_na_k! stalo
   end do

   !open (1, file='k2_n.dat')
   !do i = 1, nx
   !   write (1, '(i,2e17.8)') i, k2_res(i)
   !end do
   !close (1)
end function k2_fn

function dn_fn() result(dn_res)
   use ic, only: nx, c_double_complex, c_double, lambda, lx, im1, pi

   implicit none

   complex(c_double_complex), dimension(nx) :: dn_res
   complex(c_double_complex) k
   real(c_double) tmp
   integer i

   k = 2.0d0*pi/lambda

   dn_res(1) = dcmplx(1)

   do i = 1, nx
      tmp = 1.0d0 - (i - 1)*(i - 1)/4.0d0*(lambda/lx)*(lambda/lx)
      if (tmp >= 0) then
         dn_res(i) = dsqrt(tmp) - 1.0d0
      else
         tmp = dabs(tmp)
         dn_res(i) = -im1*dsqrt(tmp) - 1.0d0
      end if
   end do

   dn_res = k*dn_res

   open (1, file='delta_n.dat')
   do i = 1, nx
      write (1, '(i,2e17.8)') i, dn_res(i)
   end do
   close (1)
end function dn_fn

subroutine allocate_arrays()
   use ic

   implicit none

   integer(c_int) err_alloc

   allocate (a1(nx), a0(nx), ak1(nx), ak0(nx), jk1(nx), jk0(nx), atmp(nx), &
             th0(nth, 2), th1(nth, 2), dthdz(nth, 2), fk1(nx), fk2(nx), rhs0(nth, 2), z(nz), x(nx), k2(nx), &
             aktmp(nx), akzl(nx), akz0(nx), a0z0(nx), a0z0cut(nx), &
             it(it_todo), &
             ex(nx), k(2*nx), dlt(nx), tmp(nx), &
             ak_z0(nx), ak_zl(nx), a_z0(nx), a_zl(nx), & !theta(nth, 2, nz), a_amp_z0(nx), a_amp_zl(nx), sum_abs2_a_plus_by_z(nx), sum_abs2_a_plus_by_z_k(nx)&
             stat=err_alloc)

   if (err_alloc /= 0) then
      pause "allocation error"
      stop
   end if
end subroutine allocate_arrays

subroutine deallocate_arrays()
   use ic

   implicit none

   integer(c_int) err_dealloc

   deallocate (a1, a0, ak1, ak0, jk1, jk0, atmp, &
               th0, th1, dthdz, fk1, fk2, rhs0, z, x, k2, &
               aktmp, akzl, akz0, a0z0cut, &
               it, &
               ex, k, dlt, tmp, &
               !theta, a_amp_z0, a_amp_zl, sum_abs2_a_plus_by_z, sum_abs2_a_plus_by_z_k &
               stat=err_dealloc)

   if (err_dealloc /= 0) stop "deallocation error"
end subroutine deallocate_arrays

subroutine read_param() bind(c, name='read_param')
   use ic

   implicit none

   open (unit=1, file='input_fortran.in', status='old', err=101)
   read (unit=1, nml=param, err=102)
   close (unit=1)

   write (*, nml=param)

   return
101 print *, 'error of file open'; pause; stop
102 print *, 'error of reading file "input_fortran.in"'; pause; stop
end subroutine read_param

function g_fn(gx0) result(g_res)
   use ic

   implicit none

   real(c_double), dimension(nx) :: g_res
   real(c_double) gx0
   integer(c_int) i, ig0, ig1

   ig0 = gx0/hx + 1 ! dobavlyaem 1 dlia pervoy tochki
   ig1 = (nx + 1) - ig0 + 1 ! 1 k nx potomu chto furie (poslednyaya tochka ne cschitaetsya)

   if (central_mirror == .false.) then
      g_res = 1.0d0
      g_res(ig0:ig1) = 0.0d0
   else
      g_res = 0.0d0
      g_res(ig0:ig1) = 1.0d0
   end if
   end function g_fn

!subroutine step(ak0, ak1, ex, c3, jk0, jk1, dlt, h)
!   use ic, only: nx, interp, atmp
subroutine step()
   use ic

   implicit none

   !complex(8), intent(in) :: ak0(nx), ex(nx), jk0(nx), jk1(nx), dlt(nx)
   !real(8), intent(in) :: c3, h
   !complex(8), intent(out) :: ak1(nx)

   if (interp == .true.) then
      ak1(2:nx:2) = ak0(2:nx:2)*ex(2:nx:2) + &
                   c3*(jk0(2:nx:2) + jk1(2:nx:2)*(-1.0d0 + dlt(2:nx:2)*h) + &
                       ex(2:nx:2)*(jk1(2:nx:2) - jk0(2:nx:2)*(1.0d0 + dlt(2:nx:2)*h)))/dlt(2:nx:2)/dlt(2:nx:2)/h
   else
      atmp = (ak0 + c3*h/2.0d0*jk0)*cdexp(-dlt*h) !chast' a
      ak1 = atmp + c3*h/2.0d0*jk1
   end if

end subroutine step
