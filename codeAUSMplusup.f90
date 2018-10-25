program codeAUSMplusup
  !
  ! Implementation of AUSM+Up for 1D Euler equations
  ! @pcorreia
  !
  use m_init
  
  implicit none

  integer :: i

  call initialization()

  do i = 1,  niter
     call ausmplusup()
     call update()
  end do

  call output()

  print *, '** The end **'

end program codeAUSMplusup
