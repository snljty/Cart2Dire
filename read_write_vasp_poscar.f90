subroutine get_POSCAR_num_atoms(n)
    implicit none
    integer, intent(out) :: n

    character(len=6), parameter :: ifl_name = "POSCAR"
    integer :: ifl
    character(len=128) :: line
    integer :: pos
    integer :: num_atoms_this_element

    open(newunit = ifl, file = trim(ifl_name), status = "old", action = "read")
    read(ifl, '(a)') line ! comment
    read(ifl, '(a)') line ! scaling factor
    read(ifl, '(a)') line ! lattice vector 1
    read(ifl, '(a)') line ! lattice vector 2
    read(ifl, '(a)') line ! lattice vector 3
    read(ifl, '(a)') line ! may be species names?
    pos = verify(line, " ")
    if (ichar('A') <= ichar(line(pos:pos)) .and. ichar(line(pos:pos)) <= ichar('Z')) then
        read(ifl, '(a)') line ! ions per species
    end if
    close(ifl)
    n = 0
    do while (.true.)
        pos = verify(line, " ")
        if (pos == 0) then
            exit
        end if
        read(line, *) num_atoms_this_element
        pos = pos + index(line(pos:), " ") - 2
        line(:pos) = " "
        n = n + num_atoms_this_element
    end do

    return
end subroutine get_POSCAR_num_atoms

subroutine read_POSCAR(n, cell, names, coordinates, is_fraction)
    implicit none
    integer, parameter :: num_dims = 3
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(out) :: cell
    character(len=2), dimension(n), intent(out) :: names
    double precision, dimension(num_dims, n), intent(out) :: coordinates
    logical, intent(out) :: is_fraction

    external :: dscal
    double precision :: scal
    character(len=6), parameter :: ifl_name = "POSCAR", ifl_pot_name = "POTCAR"
    integer :: ifl, ifl_pot
    character(len=128) :: line, line_names
    integer :: num_atoms_this_element
    integer :: num_atoms_names_handled
    character(len=2) :: atom_name_this_element
    integer :: i
    integer :: pos

    open(newunit = ifl, file = trim(ifl_name), status = "old", action = "read")
    read(ifl, '(a)') line ! comment
    read(ifl, '(a)') line ! scaling factor
    read(line, *) scal
    read(ifl, '(a)') line ! lattice vector 1
    read(line, *) cell(:, 1)
    read(ifl, '(a)') line ! lattice vector 2
    read(line, *) cell(:, 2)
    read(ifl, '(a)') line ! lattice vector 3
    read(line, *) cell(:, 3)
    call dscal(num_dims * num_dims, scal, cell, 1)
    scal = 1.d0

    num_atoms_names_handled = 0
    read(ifl, '(a)') line_names ! may be species names?
    pos = verify(line_names, " ")
    if (ichar('A') <= ichar(line_names(pos:pos)) .and. ichar(line_names(pos:pos)) <= ichar('Z')) then
        ! contains the line of species names
        read(ifl, '(a)') line ! ions per species
        do while (.true.)
            pos = verify(line_names, " ")
            if (pos == 0) then
                exit
            end if
            read(line_names, *) atom_name_this_element
            pos = pos + index(line_names(pos:), " ") - 2
            line_names(:pos) = " "
            pos = verify(line, " ")
            read(line, *) num_atoms_this_element
            pos = pos + index(line(pos:), " ") - 2
            line(:pos) = " "
            names(num_atoms_names_handled + 1:num_atoms_names_handled + num_atoms_this_element) = atom_name_this_element
            num_atoms_names_handled = num_atoms_names_handled + num_atoms_this_element
        end do
    else
        ! does not contain the line of species names
        line = line_names(:)
        open(newunit = ifl_pot, file = trim(ifl_pot_name), status = "old", action = "read")
        do while (.true.)
            pos = verify(line, " ")
            if (pos == 0) then
                exit
            end if
            read(line, *) num_atoms_this_element
            pos = pos + index(line(pos:), " ") - 2
            line(:pos) = " "
            read(ifl_pot, '(a)') line_names
            pos = verify(line_names, " ")
            pos = pos + index(line_names(pos:), " ") - 2
            line_names(:pos) = " "
            read(line_names, *) atom_name_this_element
            do while (.true.)
                read(ifl_pot, '(a)') line_names
                if (index(line_names, "End of Dataset") /= 0) then
                    exit
                end if
            end do
            names(num_atoms_names_handled + 1:num_atoms_names_handled + num_atoms_this_element) = atom_name_this_element
            num_atoms_names_handled = num_atoms_names_handled + num_atoms_this_element
        end do
        close(ifl_pot)
    end if

    read(ifl, '(a)') line ! selective dynamics?
    pos = verify(line, " ")
    if (ichar('a') <= ichar(line(pos:pos)) .and. ichar(line(pos:pos)) <= ichar('z')) then
        line(pos:pos) = achar(ichar(line(pos:pos)) - ichar('a') + ichar('A'))
    end if
    if (line(pos:pos) /= "S") then
        backspace(ifl)
        ! does not contain selective dynamic
    end if
    read(ifl, '(a)') line ! method of coordinates
    pos = verify(line, " ")
    if (ichar('a') <= ichar(line(pos:pos)) .and. ichar(line(pos:pos)) <= ichar('z')) then
        line(pos:pos) = achar(ichar(line(pos:pos)) - ichar('a') + ichar('A'))
    end if
    if (line(pos:pos) == "C" .or. line(pos:pos) == "K") then
        is_fraction = .false.
    else
        is_fraction = .true.
    end if

    do i = 1, n
        read(ifl, '(a)') line
        read(line, *) coordinates(:, i)
    end do
    close(ifl)

    return
end subroutine read_POSCAR

subroutine write_POSCAR(n, cell, names, coordinates, is_fraction)
    implicit none
    integer, parameter :: num_dims = 3
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(in) :: cell
    character(len=2), dimension(n), intent(in) :: names
    double precision, dimension(num_dims, n), intent(in) :: coordinates
    logical, intent(in) :: is_fraction

    character(len=10), parameter :: ofl_name = "POSCAR_new"
    integer :: ofl
    integer :: i
    integer :: ind_element
    character(len=2) :: atom_name_this_element
    character(len=128) :: line_names, line
    integer :: num_atoms_this_element
    character(len=7) :: str

    open(newunit = ofl, file = trim(ofl_name), status = "replace", action = "write")
    write(ofl, '(a)') "new"
    write(ofl, '(f3.1)') 1.d0
    do i = 1, num_dims
        write(ofl, '(3(f16.8))') cell(:, i)
    end do
    atom_name_this_element = names(1)(:)
    num_atoms_this_element = 0
    line = " "
    line_names = " "
    do i = 1, n
        if (names(i) /= atom_name_this_element) then
            write(str, '(4x, i3)') num_atoms_this_element
            line_names = trim(line_names) // repeat(" ", 5) // adjustr(atom_name_this_element)
            line = trim(line) // str
            atom_name_this_element = names(i)(:)
            num_atoms_this_element = 0
        end if
        num_atoms_this_element = num_atoms_this_element + 1
    end do
    write(str, '(4x, i3)') num_atoms_this_element
    line_names = trim(line_names) // repeat(" ", 5) // adjustr(atom_name_this_element)
    line = trim(line) // str
    write(ofl, '(a)') line_names
    write(ofl, '(a)') line
    if (is_fraction) then
        write(ofl, '(a)') "Direct"
    else
        write(ofl, '(a)') "Cartesian"
    end if
    do i = 1, n
        write(ofl, '(3(f16.8))') coordinates(:, i)
    end do
    close(ofl)

    return
end subroutine write_POSCAR

subroutine convert_VASP()
    implicit none
    integer, parameter :: num_dims = 3
    integer :: n
    double precision, dimension(num_dims, num_dims) :: cell
    character(len=2), dimension(:), allocatable :: names
    double precision, dimension(:, :), allocatable :: coordinates, coordinates_new
    logical :: is_fraction

    external :: get_POSCAR_num_atoms, read_POSCAR, write_POSCAR
    external :: Dire2Cart, Cart2Dire

    call get_POSCAR_num_atoms(n)
    allocate(names(n))
    allocate(coordinates(num_dims, n), coordinates_new(num_dims, n))
    call read_POSCAR(n, cell, names, coordinates, is_fraction)
    if (is_fraction) then
        call Dire2Cart(n, cell, coordinates, coordinates_new)
    else
        call Cart2Dire(n, cell, coordinates, coordinates_new)
    end if
    call write_POSCAR(n, cell, names, coordinates_new, .not. is_fraction)
    deallocate(names)
    deallocate(coordinates, coordinates_new)

    return
end subroutine convert_VASP

