subroutine locate_label(line, label, pos_equal)
    implicit none
    character(len=*), intent(in) :: line, label
    integer, intent(out) :: pos_equal

    integer :: pos1, pos2

    pos_equal = 0

    pos1 = index(line, label)
    if (pos1 == 0) then
        return
    end if
    pos2 = verify(line, " ")
    if (pos2 == 0) then
        return
    end if
    if (pos1 /= pos2) then
        return
    end if
    pos_equal = index(line, "=")
    if (pos_equal == 0) then
        return
    end if
    if (line(pos1 + len_trim(label):pos_equal - 1) == " ") then
        return
    end if
    pos_equal = 0
    return
end subroutine locate_label

subroutine get_pw_inp_num_atoms(ifl_name, n)
    implicit none
    character(len=*), intent(in) :: ifl_name
    integer, intent(out) :: n
    integer :: ifl
    character(len=128) :: line
    integer :: pos

    external :: locate_label

    ! actually all these should be in the "&system" block.

    open(newunit = ifl, file = trim(ifl_name), status = "old", action = "read")
    do while (.true.)
        read(ifl, '(a)') line
        call locate_label(line, "nat", pos)
        if (pos /= 0) then
            exit
        end if
    end do
    close(ifl)
    read(line(pos + 1:), *) n

    return
end subroutine get_pw_inp_num_atoms

subroutine read_pw_inp(ifl_name, n, cell, names, coordinates, is_fraction)
    implicit none
    integer, parameter :: num_dims = 3
    character(len=*), intent(in) :: ifl_name
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(out) :: cell
    character(len=2), dimension(n), intent(out) :: names
    double precision, dimension(num_dims, n), intent(out) :: coordinates
    logical, intent(out) :: is_fraction

    external :: locate_label
    external :: dscal

    double precision, parameter :: Bohr2Angstrom = 5.2917721092d-1
    integer :: ifl
    integer :: i, j
    character(len=128) :: line
    integer :: pos
    character(len=16) :: unit_name
    integer :: status
    integer :: ibrav
    logical :: is_celldm_provided, is_cell_para_provided
    double precision :: celldm(2 * num_dims), A, B, C, cosAB, cosAC, cosBC
    double precision :: alat
    character(len=9) :: str
    double precision :: c_div_a, b_div_a, cos_ang, tmp(num_dims)

    ! actually all these should be in the "&system" block.

    is_celldm_provided = .false.
    is_cell_para_provided = .false.
    celldm = 0.d0
    A = 0.d0
    B = 0.d0
    C = 0.d0
    cosAB = 0.d0
    cosAC = 0.d0
    cosBC = 0.d0
    open(newunit = ifl, file = trim(ifl_name), status = "old", action = "read")
    ! get ibrav
    do while (.true.)
        read(ifl, '(a)', iostat = status) line
        if (status /= 0) then
            write(0, '(a)') 'Error! Lack of "ibrav" is currently unsupported.'
            close(ifl)
            stop 1
        end if
        call locate_label(line, "ibrav", pos)
        if (pos /= 0) then
            exit
        end if
    end do
    read(line(pos + 1:), *) ibrav

    ! get celldm if possible
    rewind(ifl)
    do while (.true.)
        read(ifl, '(a)', iostat = status) line
        if (status /= 0) then
            exit
        end if
        i = 1
        write(str, '(a, i0, a)') "celldm(", i, ")"
        call locate_label(line, str, pos)
        if (pos /= 0) then
            is_celldm_provided = .true.
            read(line(pos + 1:), *) celldm(1)
            do i = 2, 2 * num_dims
                rewind(ifl)
                do while (.true.)
                    read(ifl, '(a)', iostat = status) line
                    if (status /= 0) then
                        exit
                    end if
                    write(str, '(a, i0, a)') "celldm(", i, ")"
                    call locate_label(line, trim(str), pos)
                    if (pos /= 0) then
                        read(line(pos + 1:), *) celldm(i)
                        exit
                    end if
                end do
            end do
            exit
        end if
    end do

    ! get A B C cosAB cosAC cosBC if possible
    rewind(ifl)
    do while (.true.)
        read(ifl, '(a)', iostat = status) line
        if (status /= 0) then
            exit
        end if
        call locate_label(line, "A", pos)
        if (pos /= 0) then
            is_cell_para_provided = .true.
            read(line(pos + 1:), *) A
            exit
        end if
    end do
    if (is_cell_para_provided) then
        rewind(ifl)
        do while (.true.)
            read(ifl, '(a)', iostat = status) line
            if (status /= 0) then
                exit
            end if
            call locate_label(line, "B", pos)
            if (pos /= 0) then
                read(line(pos + 1:), *) B
                exit
            end if
        end do
        rewind(ifl)
        do while (.true.)
            read(ifl, '(a)', iostat = status) line
            if (status /= 0) then
                exit
            end if
            call locate_label(line, "C", pos)
            if (pos /= 0) then
                read(line(pos + 1:), *) C
                exit
            end if
        end do
        rewind(ifl)
        do while (.true.)
            read(ifl, '(a)', iostat = status) line
            if (status /= 0) then
                exit
            end if
            call locate_label(line, "cosAB", pos)
            if (pos /= 0) then
                read(line(pos + 1:), *) cosAB
                exit
            end if
        end do
        rewind(ifl)
        do while (.true.)
            read(ifl, '(a)', iostat = status) line
            if (status /= 0) then
                exit
            end if
            call locate_label(line, "cosAC", pos)
            if (pos /= 0) then
                read(line(pos + 1:), *) cosAC
                exit
            end if
        end do
        rewind(ifl)
        do while (.true.)
            read(ifl, '(a)', iostat = status) line
            if (status /= 0) then
                exit
            end if
            call locate_label(line, "cosBC", pos)
            if (pos /= 0) then
                read(line(pos + 1:), *) cosBC
                exit
            end if
        end do
    end if

    ! cannot specify celldm and A B C cosAB cosAC cosBC together.
    if (is_celldm_provided .and. is_cell_para_provided) then
        write(0, '(a)') 'Error! Cannot provide "celldm(1)" and "A" together.'
        close(ifl)
        stop 1
    end if
    if (is_celldm_provided) then
        alat = celldm(1) * Bohr2Angstrom
    end if
    if (is_cell_para_provided) then
        alat = A
    end if

    ! read cell
    cell = 0.d0
    rewind(ifl)
    ! read from CELL_PARAMETERS
    if (ibrav == 0) then
        do while (.true.)
            read(ifl, '(a)') line
            pos = index(line, "CELL_PARAMETERS")
            if (pos /= 0) then
                exit
            end if
        end do
        line(pos:pos + len("CELL_PARAMETERS") - 1) = " "
        read(line, *, iostat = status) unit_name
        if (status /= 0) then
            unit_name = " "
        end if
        if (unit_name == " ") then
            write(0, "(a)") "Warning: non-specified lattice parameter unit is deprecated."
        end if
        do i = 1, 3
            read(ifl, '(a)') line
            read(line, *) cell(:, i)
        end do
        if (unit_name == " ") then
            if (is_celldm_provided .or. is_cell_para_provided) then
                unit_name = "alat"
            else
                unit_name = "bohr"
            end if
        end if
        if (unit_name == "alat") then
            if (is_celldm_provided .or. is_cell_para_provided) then
                call dscal(num_dims * num_dims, alat, cell, 1)
                alat = 1.d0
            else
                write(0, '(a)') "Error! In this situation either celldm(1) or A should be specified"
                close(ifl)
                stop 1
            end if
        else if (unit_name == "bohr") then
            call dscal(num_dims * num_dims, Bohr2Angstrom, cell, 1)
        else if (unit_name == "angstrom") then
            continue
        else
            write(0, '(a)') 'Error! Unrecognized cell unit: "', trim(unit_name), '".'
            close(ifl)
            stop 1
        end if
    ! Cubic P (sc)
    else if (ibrav == 1) then
        cell(1, 1) = 1.d0
        cell(2, 2) = 1.d0
        cell(3, 3) = 1.d0
    ! Cubic F (fcc)
    else if (ibrav == 2) then
        cell(1, 1) = - 5.d-1
        cell(3, 1) = 5.d-1
        cell(2, 2) = 5.d-1
        cell(3, 2) = 5.d-1
        cell(1, 3) = - 5.d-1
        cell(2, 3) = 5.d-1
    ! Cubic I (bcc)
    else if (ibrav == 3) then
        cell = 5.d-1
        cell(1, 2) = - 5.d-1
        cell(1, 3) = - 5.d-1
        cell(2, 3) = - 5.d-1
    ! Cubic I (bcc), more symmetric axis
    else if (ibrav == -3) then
        cell = 5.d-1
        cell(1, 1) = - 5.d-1
        cell(2, 2) = - 5.d-1
        cell(3, 3) = - 5.d-1
    ! Hexagonal and Trigonal P
    else if (ibrav == 4) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
        else
            c_div_a = C / A
        end if
        cell(1, 1) = 1.d0
        cell(1, 2) = -5.d-1
        cell(2, 2) = sqrt(3.d0) / 2.d0
        cell(3, 3) = c_div_a
    ! Trigonal R, 3fold axis c
    else if (ibrav == 5) then
        if (is_celldm_provided) then
            cos_ang = celldm(4)
        else
            cos_ang = cosAB ! cosBC?
        end if
        tmp(1)=sqrt((1.d0 - cos_ang) / 2.d0)
        tmp(2)=sqrt((1.d0 - cos_ang) / 6.d0)
        tmp(3)=sqrt((1.d0 + 2.d0 * cos_ang) / 3.d0)
        cell(1, 1) = tmp(1)
        cell(2, 1) = - tmp(2)
        cell(3, 1) = tmp(3)
        cell(2, 2) = 2 * tmp(2)
        cell(3, 2) = tmp(3)
        cell(1, 3) = - tmp(1)
        cell(2, 3) = - tmp(2)
        cell(3, 3) = tmp(3)
    ! Trigonal R, 3fold axis <111>
    else if (ibrav == -5) then
        if (is_celldm_provided) then
            cos_ang = celldm(4)
        else
            cos_ang = cosAB ! cosBC?
        end if
        tmp(1)=sqrt((1.d0 - cos_ang) / 2.d0)
        tmp(2)=sqrt((1.d0 - cos_ang) / 6.d0)
        tmp(3)=sqrt((1.d0 + 2.d0 * cos_ang) / 3.d0)
        b_div_a = tmp(3) - 2.d0 * sqrt(2.d0) * tmp(2)
        c_div_a = tmp(3) + sqrt(2.d0) * tmp(2)
        do i = 1, num_dims
            do j = 1, num_dims
                if (i == j) then
                    cell(i, j) = tmp(3) - 2.d0 * sqrt(2.d0) * tmp(2)
                    ! cell(i, j) = tmp(3) + 2.d0 * sqrt(2.d0) * tmp(2)
                else
                    cell(i, j) = tmp(3) + sqrt(2.d0) * tmp(2)
                    ! cell(i, j) = tmp(3) - sqrt(2.d0) * tmp(2)
                end if
            end do
        end do
    ! Tetragonal P (st)
    else if (ibrav == 6) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
        else
            c_div_a = C / A
        end if
        cell(1, 1) = 1.d0
        cell(2, 2) = 1.d0
        cell(3, 3) = c_div_a
    ! Tetragonal I (bct)
    else if (ibrav == 7) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
        else
            c_div_a = C / A
        end if
        cell(1, 1) = 5.d-1
        cell(2, 1) = -5.d-1
        cell(3, 1) = c_div_a
        cell(1, 2) = 5.d-1
        cell(2, 2) = 5.d-1
        cell(3, 2) = c_div_a
        cell(1, 3) = -5.d-1
        cell(2, 3) = -5.d-1
        cell(3, 3) = c_div_a
    ! Orthorhombic P
    else if (ibrav == 8) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 1.d0
        cell(2, 2) = b_div_a
        cell(3, 3) = c_div_a
    ! Orthorhombic base-centered(bco)
    else if (ibrav == 9) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 5.d-1
        cell(2, 1) = b_div_a / 2.d0
        cell(1, 2) = -5.d-1
        cell(2, 2) = b_div_a / 2.d0
        cell(3, 3) = c_div_a
    ! as 9, alternate description
    else if (ibrav == -9) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 5.d-1
        cell(2, 1) = - b_div_a / 2.d0
        cell(1, 2) = 5.d-1
        cell(2, 2) = b_div_a / 2.d0
        cell(3, 3) = c_div_a
    ! Orthorhombic one-face base-centered A-type
    else if (ibrav == 91) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(2, 2) = b_div_a / 2.d0
        cell(3, 2) = - c_div_a / 2.d0
        cell(2, 3) = b_div_a / 2.d0
        cell(3, 3) = c_div_a / 2.d0
    ! Orthorhombic face-centered
    else if (ibrav == 10) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 5.d-1
        cell(3, 1) = c_div_a / 2.d0
        cell(1, 2) = 5.d-1
        cell(2, 2) = b_div_a / 2.d0
        cell(2, 3) = b_div_a / 2.d0
        cell(3, 3) = c_div_a / 2.d0
    ! Orthorhombic body-centered
    else if (ibrav == 11) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 5.d-1
        cell(2, 1) = b_div_a / 2.d0
        cell(3, 1) = c_div_a / 2.d0
        cell(1, 2) = -5.d-1
        cell(2, 2) = b_div_a / 2.d0
        cell(3, 2) = c_div_a / 2.d0
        cell(1, 3) = -5.d-1
        cell(2, 3) = - b_div_a / 2.d0
        cell(3, 3) = - c_div_a / 2.d0
    ! Monoclinic P, unique axis c
    else if (ibrav == 12) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
            cos_ang = celldm(4) ! a and b
        else
            c_div_a = C / A
            b_div_a = B / A
            cos_ang = cosAB
        end if
        cell(1, 1) = 1.d0
        cell(1, 2) = b_div_a * cos_ang
        cell(2, 2) = b_div_a * sqrt(1.d0 - cos_ang ** 2)
        cell(3, 3) = c_div_a
    ! Monoclinic P, unique axis b
    else if (ibrav == -12) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
            cos_ang = celldm(5) ! a and c
        else
            c_div_a = C / A
            b_div_a = B / A
            cos_ang = cosAC
        end if
        cell(1, 1) = 1.d0
        cell(2, 2) = b_div_a
        cell(1, 3) = c_div_a * cos_ang
        cell(3, 3) = c_div_a * sqrt(1.d0 - cos_ang ** 2)
    ! Monoclinic base-centered (unique axis c)
    else if (ibrav == 13) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
            cos_ang = celldm(4) ! a and b
        else
            c_div_a = C / A
            b_div_a = B / A
            cos_ang = cosAB
        end if
        cell(1, 1) = 5.d-1
        cell(3, 1) = - c_div_a / 2.d0
        cell(1, 2) = b_div_a * cos_ang
        cell(2, 2) = b_div_a * sqrt(1.d0 - cos_ang ** 2)
        cell(1, 3) = 5.d-1
        cell(3, 3) = c_div_a
    ! Monoclinic base-centered (unique axis b)
    else if (ibrav == -13) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
            cos_ang = celldm(5) ! a and c
        else
            c_div_a = C / A
            b_div_a = B / A
            cos_ang = cosAC
        end if
        cell(1, 1) = 5.d-1
        cell(2, 1) = b_div_a / 2.d0
        cell(1, 2) = -5.d-1
        cell(2, 2) = b_div_a / 2.d0
        cell(1, 3) = c_div_a * cos_ang
        cell(3, 3) = c_div_a * sqrt(1.d0 - cos_ang ** 2)
    ! Triclinic
    else if (ibrav == 14) then
        if (is_celldm_provided) then
            c_div_a = celldm(3)
            b_div_a = celldm(2)
        else
            c_div_a = C / A
            b_div_a = B / A
        end if
        cell(1, 1) = 1.d0
        if (is_celldm_provided) then
            cell(1, 2) = b_div_a * celldm(6)
            cell(2, 2) = b_div_a * sqrt(1.d0 - celldm(6) ** 2)
            cell(1, 3) = c_div_a * celldm(5)
            cell(2, 3) = c_div_a * (celldm(4) - celldm(5) * celldm(6)) / sqrt(1.d0 - celldm(6) ** 2)
            cell(3, 3) = c_div_a * sqrt(1.d0 + 2.d0 * celldm(4) * celldm(5) * celldm(6) &
                    - celldm(4) ** 2 - celldm(5) ** 2 - celldm(6) ** 2) / sqrt(1.d0 - celldm(6) ** 2)
        else
            cell(1, 2) = b_div_a * cosAB
            cell(2, 2) = b_div_a * sqrt(1.d0 - cosAB ** 2)
            cell(1, 3) = c_div_a * cosAC
            cell(2, 3) = c_div_a * (cosBC - cosAC * cosAB) / sqrt(1.d0 - cosAB ** 2)
            cell(3, 3) = c_div_a * sqrt(1.d0 + 2.d0 * cosBC * cosAC * cosAB &
                    - cosBC ** 2 - cosAC ** 2 - cosAB ** 2) / sqrt(1.d0 - cosAB ** 2)
        end if
    else
        write(0, '(a, i0, a)') "Error! Unrecognized ibrav = ", ibrav, "."
        close(ifl)
        stop 1
    end if
    if (ibrav /= 0) then
        call dscal(num_dims * num_dims, alat, cell, 1)
        alat = 1.d0
    end if

    ! read atomic names and coordinates
    rewind(ifl)
    do while (.true.)
        read(ifl, '(a)') line
        pos = index(line, "ATOMIC_POSITIONS")
        if (pos /= 0) then
            exit
        end if
    end do
    line(pos:pos + len("ATOMIC_POSITIONS") - 1) = " "
    read(line, *, iostat = status) unit_name
    if (status /= 0) then
        unit_name = " "
    end if
    if (unit_name == " ") then
        unit_name = "alat"
        write(0, "(a)") "Warning: non-specified atomic positions units are deprecated."
    end if
    if (unit_name == "crystal_sg") then
        write(0, '(a)') 'Error: atomic positions units "crystal_sg" is currently unsupported.'
        close(ifl)
        stop 1
    end if
    do i = 1, n
        read(ifl, '(a)') line
        read(line, *) names(i), coordinates(:, i)
    end do
    close(ifl)
    if (unit_name == "crystal") then
        is_fraction = .true.
    else if (unit_name == "bohr") then
        is_fraction = .false.
        call dscal(n * num_dims, Bohr2Angstrom, coordinates, 1)
    else if (unit_name == "alat") then
        is_fraction = .false.
        call dscal(n * num_dims, alat, coordinates, 1)
    else if (unit_name == "angstrom") then
        is_fraction = .false.
        continue
    else
        write(0, '(a, a, a)') 'Error! Unrecognized atomic positions units: "', trim(unit_name), '".'
        stop 1
    end if

    return
end subroutine read_pw_inp

subroutine write_pw_inp(ofl_name, n, cell, names, coordinates, is_fraction)
    implicit none
    integer, parameter :: num_dims = 3
    character(len=*), intent(in) :: ofl_name
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(in) :: cell
    character(len=2), dimension(n), intent(in) :: names
    double precision, dimension(num_dims, n), intent(in) :: coordinates
    logical, intent(in) :: is_fraction

    integer :: ofl
    integer :: i

    open(newunit = ofl, file = trim(ofl_name), status = "replace", action = "write")
    write(ofl, '(a)') "&system"
    write(ofl, '(1x, a, i0)') "ibrav = ", 0
    write(ofl, '(1x, a, i0)') "nat = ", n
    write(ofl, '(a)') "/"

    write(ofl, '(a)') "CELL_PARAMETERS angstrom"
    do i = 1, num_dims
        write(ofl, '(3(f16.8))') cell(:, i)
    end do

    write(ofl, '(a)', advance = "no") "ATOMIC_POSITIONS "
    if (is_fraction) then
        write(ofl, '(a)') "crystal"
    else
        write(ofl, '(a)') "angstrom"
    end if
    do i = 1, n
        write(ofl, '(1x, a2, 3(f16.8))') names(i), coordinates(:, i)
    end do
    close(ofl)

    return
end subroutine write_pw_inp

subroutine convert_QE(ifl_name)
    implicit none
    character(len=*), intent(in) :: ifl_name

    character(len=128) :: ofl_name
    integer, parameter :: num_dims = 3
    integer :: n
    double precision, dimension(num_dims, num_dims) :: cell
    character(len=2), dimension(:), allocatable :: names
    double precision, dimension(:, :), allocatable :: coordinates, coordinates_new
    logical :: is_fraction
    integer :: pos

    external :: get_POSCAR_num_atoms, read_POSCAR, write_POSCAR
    external :: Dire2Cart, Cart2Dire

    pos = index(ifl_name, ".", back = .true.)
    if (pos == 0) then
        ofl_name = trim(ifl_name) // "_new"
    else
        ofl_name = ifl_name(:pos - 1) // "_new" // trim(ifl_name(pos:))
    end if

    call get_pw_inp_num_atoms(ifl_name, n)
    allocate(names(n))
    allocate(coordinates(num_dims, n), coordinates_new(num_dims, n))
    call read_pw_inp(ifl_name, n, cell, names, coordinates, is_fraction)
    if (is_fraction) then
        call Dire2Cart(n, cell, coordinates, coordinates_new)
    else
        call Cart2Dire(n, cell, coordinates, coordinates_new)
    end if
    call write_pw_inp(ofl_name, n, cell, names, coordinates_new, .not. is_fraction)
    deallocate(names)
    deallocate(coordinates, coordinates_new)

    return
end subroutine convert_QE

