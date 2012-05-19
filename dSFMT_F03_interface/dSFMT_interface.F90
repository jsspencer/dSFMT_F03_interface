module dSFMT_interface

! copyright James Spencer 2012.
! New BSD License, see License.txt for details.

! This contains a wrapper around a subset of of the functionality in dSFMT RNG
! (double precision SIMD-oriented Fast Mersenne Twister random number
! generator).

! Functions which return a single (pseudo-)random number are defined as inline
! in the original dSFMT source.  Accessing these functions is quite slow
! (compared to C) in Fortran as we can't inline them (or indeed call them
! without modifying the C source).  A more efficient way is to fill up an array
! of random numbers and then access the array as a stream, refilling the array
! as needed.  By thus amortizing the cost of calling non-inlined functions
! across many random numbers, the cost becomes very competitive with
! a C implementation.  Functions are provided which return a single random
! number and refill the storage array when required, thus giving a fast
! interface which closely mirrors that of the C implementation.

! Usage (for each psuedo-random stream):
!   i)   call dSFMT_init.
!   ii)  access desired number of random numbers either by calling
!        a get_rand_xxx function the required number of times or by filling up
!        user-defined arrays with a dsfmt_fill_array_xxx subroutine.
!   iii) call dSFMT_end to free all memory.

! IMPORTANT: see the warning below about calling multiple get_rand_xxx
! functions.

! Match definition of DSFMT_MEXP in dSFMT.h.
#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
#warning "DSFMT_MEXP is not defined. I assume DSFMT_MEXP is 19937."
#endif
#define DSFMT_MEXP 19937
#endif

use, intrinsic :: iso_c_binding
implicit none

private
public :: dSFMT_t, dSFMT_init, dSFMT_end, dSFMT_reset,  &
          fill_array_close_open, fill_array_open_close, &
          fill_array_open_open, get_rand_close_open,    &
          get_rand_open_close, get_rand_open_open

! Expose functions from C as needed.
! See dSFMT documentation for details.
interface
    pure function dsfmt_get_min_array_size() result(min_size) bind(c, name='dsfmt_get_min_array_size')
        import :: c_int32_t
        integer(c_int32_t) :: min_size
    end function dsfmt_get_min_array_size
    pure function malloc_dsfmt_t() result(dsfmt_state) bind(c, name='malloc_dsfmt_t')
        import :: c_ptr
        type(c_ptr) :: dsfmt_state
    end function malloc_dsfmt_t
    pure subroutine free_dsfmt_t(dsfmt_state) bind(c, name='free_dsfmt_t')
        import :: c_ptr
        type(c_ptr), value, intent(in) :: dsfmt_state
    end subroutine free_dsfmt_t
    pure subroutine dsfmt_chk_init_gen_rand(dsfmt_state, seed, mexp) bind(c, name='dsfmt_chk_init_gen_rand')
        import :: c_ptr, c_int32_t
        type(c_ptr), value, intent(in) :: dsfmt_state
        integer(c_int32_t), value, intent(in) :: seed, mexp
    end subroutine dsfmt_chk_init_gen_rand
    pure subroutine dsfmt_fill_array_close_open(dSFMT_state, array, array_size) bind(c, name='dsfmt_fill_array_close_open')
        import :: c_ptr, c_double, c_int32_t
        type(c_ptr), value, intent(in) :: dSFMT_state
        real(c_double), intent(out) :: array(*)
        integer(c_int32_t), value, intent(in) :: array_size
    end subroutine dsfmt_fill_array_close_open
    pure subroutine dsfmt_fill_array_open_close(dSFMT_state, array, array_size) bind(c, name='dsfmt_fill_array_open_close')
        import :: c_ptr, c_double, c_int32_t
        type(c_ptr), value, intent(in) :: dSFMT_state
        real(c_double), intent(out) :: array(*)
        integer(c_int32_t), value, intent(in) :: array_size
    end subroutine dsfmt_fill_array_open_close
    pure subroutine dsfmt_fill_array_open_open(dSFMT_state, array, array_size) bind(c, name='dsfmt_fill_array_open_open')
        import :: c_ptr, c_double, c_int32_t
        type(c_ptr), value, intent(in) :: dSFMT_state
        real(c_double), intent(out) :: array(*)
        integer(c_int32_t), value, intent(in) :: array_size
    end subroutine dsfmt_fill_array_open_open
end interface

type dSFMT_t
    private
    ! Pointer to the dsfmt_t internal state (as defined in dSFMT.h).
    type(c_ptr) :: dSFMT_state
    ! Testing indicates that 50000 is a very good size for the array storing the
    ! random numbers.  Testing was done standalone, so undoubtedly influenced by
    ! cache size and this might be different for real-world applications, but it's easy to
    ! change to allocatable later on.
    integer(c_int) :: random_store_size
    ! Seed passed to dSFMT.
    integer :: seed
    ! Store of random numbers.
    real(c_double), allocatable :: random_store(:)
    ! Index of the next random number to be returned from random_store.
    integer :: next_element
end type dSFMT_t

integer, parameter, private :: dp = selected_real_kind(15,307)

contains

!--- Initialisation, termination and state interaction ---

    pure subroutine dSFMT_init(seed, rng_store_size, rng)

        ! Initialise the dSFMT RNG and fill rng%random_store with
        ! a block of random numbers in interval [0,1).
        !
        ! In:
        !    seed: seed for the RNG.
        !    rng_store_size: number of random numbers to store at once.
        ! Out:
        !    rng: dSFMT_t with internal variables initialised and associated
        !       with an initialised psuedo-random number stream.

        integer, intent(in) :: seed, rng_store_size
        type(dSFMT_t), intent(inout) :: rng

        if (dsfmt_get_min_array_size() > rng_store_size) then
            rng%random_store_size = dsfmt_get_min_array_size()
        else
            rng%random_store_size = rng_store_size
        end if

        allocate(rng%random_store(rng%random_store_size))

        ! Set current element to be larger than the store, so it is
        ! filled in the first call to the get_* functions.
        rng%next_element = rng%random_store_size + 1

        rng%seed = int(seed, c_int32_t)

        ! Create dSFMT state
        rng%dSFMT_state = malloc_dsfmt_t()

        ! Initialise dSFMT PRNG.
        call dsfmt_chk_init_gen_rand(rng%dSFMT_state, rng%seed, DSFMT_MEXP)

    end subroutine dSFMT_init

    pure subroutine dSFMT_end(rng)

        ! Deallocate and reset a dSFMT_t variable safely.

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is closed on output and memory
        !       associated with it is deallocated.

        type(dSFMT_t), intent(inout) :: rng

        ! Destroy dSFMT state
        call free_dsfmt_t(rng%dSFMT_state)

        ! Reset store.
        deallocate(rng%random_store)
        rng%random_store_size = -1
        rng%next_element = -1

    end subroutine dSFMT_end

    pure subroutine dSFMT_reset(rng)

        ! Reset the dSFMT state such that the next get_rand_xxx call requires
        ! the store of random numbers to be refilled.

        ! In/Out:
        !    rng: dSFMT_t.  The internal store of random numbers is cleared.

        type(dSFMT_t), intent(inout) :: rng

        rng%next_element = rng%random_store_size + 1
        rng%random_store = huge(1.0_dp)

    end subroutine dSFMT_reset

!--- Get arrays of random numbers ---

    pure subroutine fill_array_close_open(rng, array)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.
        ! Out:
        !    array: array filled with random numbers in the interval [0,1).

        ! NOTE: array must be at least as large as the size returned by dsfmt_get_min_array_size.

        type(dSFMT_t), intent(inout) :: rng
        real(c_double), intent(out) :: array(:)

        call dsfmt_fill_array_close_open(rng%dSFMT_state, array, size(array))

    end subroutine fill_array_close_open

    pure subroutine fill_array_open_close(rng, array)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.
        ! Out:
        !    array: array filled with random numbers in the interval (0,1].

        ! NOTE: array must be at least as large as the size returned by dsfmt_get_min_array_size.

        type(dSFMT_t), intent(inout) :: rng
        real(c_double), intent(out) :: array(:)

        call dsfmt_fill_array_open_close(rng%dSFMT_state, array, size(array))

    end subroutine fill_array_open_close

    pure subroutine fill_array_open_open(rng, array)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.
        ! Out:
        !    array: array filled with random numbers in the interval (0,1).

        ! NOTE: array must be at least as large as the size returned by dsfmt_get_min_array_size.

        type(dSFMT_t), intent(inout) :: rng
        real(c_double), intent(out) :: array(:)

        call dsfmt_fill_array_open_open(rng%dSFMT_state, array, size(array))

    end subroutine fill_array_open_open

!--- Get a random number ---

! WARNING: these calls use a store internal to the dSFMT_t variable.  Calls to
! different functions with the same state must be separated either by a call to
! dSFMT_reset or a call to dSFMT_end followed by reinitialisation of the dSFMT_t
! state using dSFMT_init.

    function get_rand_close_open(rng) result(r)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.  The store of random numbers is refilled if necessary.

        ! Returns:
        !    random number in interval [0,1).

        real(dp) :: r
        type(dSFMT_t), intent(inout) :: rng

        if (rng%next_element == rng%random_store_size+1) then
            ! Run out of random numbers: get more.
            call dsfmt_fill_array_close_open(rng%dSFMT_state, rng%random_store, rng%random_store_size)
            rng%next_element = 1
        end if

        r = rng%random_store(rng%next_element)
        rng%next_element = rng%next_element + 1 

    end function get_rand_close_open

    function get_rand_open_close(rng) result(r)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.  The store of random numbers is refilled if necessary.

        ! Returns:
        !    random number in interval (0,1].

        real(dp) :: r
        type(dSFMT_t), intent(inout) :: rng

        if (rng%next_element == rng%random_store_size+1) then
            ! Run out of random numbers: get more.
            call dsfmt_fill_array_open_close(rng%dSFMT_state, rng%random_store, rng%random_store_size)
            rng%next_element = 1
        end if

        r = rng%random_store(rng%next_element)
        rng%next_element = rng%next_element + 1 

    end function get_rand_open_close

    function get_rand_open_open(rng) result(r)

        ! In/Out:
        !    rng: dSFMT_t.  The dSFMT state is updated by the request for random
        !       numbers.  The store of random numbers is refilled if necessary.

        ! Returns:
        !    random number in interval (0,1).

        real(dp) :: r
        type(dSFMT_t), intent(inout) :: rng

        if (rng%next_element == rng%random_store_size+1) then
            ! Run out of random numbers: get more.
            call dsfmt_fill_array_open_open(rng%dSFMT_state, rng%random_store, rng%random_store_size)
            rng%next_element = 1
        end if

        r = rng%random_store(rng%next_element)
        rng%next_element = rng%next_element + 1 

    end function get_rand_open_open

end module dSFMT_interface
