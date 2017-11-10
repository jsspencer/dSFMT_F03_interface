program test

    use dSFMT_interface
    implicit none

    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: N = 10000000

    type(dSFMT_t) :: rng

    real(dp) :: rand_sum
    integer :: i
    character(:), allocatable :: state_str

    call dSFMT_init(7, 50000, rng)

    rand_sum = 0.0_dp
    do i = 1, N
        rand_sum = rand_sum + get_rand_close_open(rng)
    end do

    write (6,'(a7,i8,a19,f14.6)') 'Sum of ', N, ' random numbers is ', rand_sum

    call dsfmt_state_to_string(rng, state_str)

    do i = 0, 9
        write (6,'(i1,1x,f8.6)') i, get_rand_close_open(rng)
    end do

    call dSFMT_end(rng)

    call dSFMT_init(7, 50000, rng)
    call dsfmt_str_to_state(rng, state_str)

    do i = 0, 9
        write (6,'(i1,1x,f8.6)') i, get_rand_close_open(rng)
    end do

    call dSFMT_end(rng)
    deallocate(state_str)

end program test
