program main
    implicit none

    integer :: argc, iarg, argl_max, argl
    character(len=:), dimension(:), allocatable :: argv
    character(len=:), allocatable :: argv0

    character(len=128) :: ifl_name
    logical :: alive

    external :: convert_VASP, convert_QE

    argc = command_argument_count()
    call get_command_argument(0, length = argl)
    allocate(character(len=argl) :: argv0)
    call get_command_argument(0, argv0)
    argl_max = 0
    do iarg = 1, argc
        call get_command_argument(iarg, length = argl)
        if (argl > argl_max) then
            argl_max = argl
        end if
    end do
    allocate(character(len=argl_max) :: argv(argc))
    do iarg = 1, argc
        call get_command_argument(iarg, argv(iarg))
    end do

    if (argc >= 2) then
        write(0, '(a, a)') trim(argv0), ":"
        write(0, '(a)') "Convert between fractional crystal coordinates and Cartesian crystal coordinates."
        write(0, '(a)') "Supports input file of VASP (POSCAR) and input file of pw.x of Quantum ESPRESSO."
        write(0, '()')
        write(0, '(a, a, a)') "Usage: ", trim(argv0), " [QE_pw_input_file]"
        write(0, '()')
        write(0, '(a)') "If QE_pw_input_file is not provided, "
        write(0, '(a)') "it will try to read the POSCAR of VASP in the current directory, "
        write(0, '(a)') "and if POSCAR is not found, it will ask for the name of the Quantum ESPRESSO "
        write(0, '(a)') "pw.x input file."
        write(0, '(a)') "The output name has the same suffix (if any) with the input file, and the file name "
        write(0, '(a)') "is appended by ""_new""."
        write(0, '(a)') "if the input file is fractional crystal coordinates (a.k.a. direct coordinates), "
        write(0, '(a)') "the output file will be Cartesian coordinates, and vice versa."
        write(0, '()')
        write(0, '(a)') "Exiting now."
        stop 1
    else if (argc == 1) then
        ifl_name = argv(1)(:)
        call convert_QE(ifl_name)
    else
        inquire(file = "POSCAR", exist = alive)
        if (alive) then
            call convert_VASP()
        else
            write(*, '(a)') "Input file name of input file of pw.x of Quantum ESPRESSO:"
            read(*, '(a)') ifl_name
            if (ifl_name(1:1) == '"') then
                ifl_name = ifl_name(2:)
                if (ifl_name(len_trim(ifl_name):) == '"') then
                    ifl_name(len_trim(ifl_name):) = " "
                end if
            end if
            call convert_QE(ifl_name)
        end if
    end if

    write(*, '(a)') "Done."

    deallocate(argv0)
    deallocate(argv)

    stop
end program main

