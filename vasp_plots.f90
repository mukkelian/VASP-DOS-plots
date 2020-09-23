	implicit none
        real, allocatable :: x(:), y(:), z(:) !x,y,z-co-ordinates of ions
        integer, allocatable :: ion(:)
	real, allocatable :: a(:,:,:), b(:,:,:), hdr(:,:), tot(:,:)
	real wtt
        character*20 filename,coordinate,title,car, sp_name, ttl
        character(len=2),allocatable :: species(:)

        character(len=11), dimension(70) :: str
        character(len=11), dimension(1) :: str1
	character orb, so
	integer orb_l, ol, jj
        integer num,line,line1,i,j,n_atoms, n_species, ispin, k
        integer at1, at2
        real a1,b1,c1,d1,a2,b2,c2,d2,e2,a3, abc(3,3)
        real  emax, emin, eferm, wt
        num = 0
	!_______________ INPUT FILE __________________!
        open(unit=0,file='input',status='old',action='read')
	read(0,*) n_species ! no. of species
	read(0,*) ispin ! type of spin calc
	read(0,*) orb	! last family of orbitals in PROCAR, i.e., s,p,d or f
	read(0,*) so	! Spin-Orbit coupling
	close (0)

	print*,'!'
	print*,'!'
	print*, '~ Mukesh Kumar Sharma'
	print*, 'contact:: msharma1@ph.iitr.ac.in'
	print*,'!'
	print*,'!#########################!'
	print*, '!May the Force be with you!'
	print*,'!#########################!'
	print*,'.'

	open(10001,file='POSCAR',status='unknown')

	!____________ READING STRUCTURE FILE _________!
	read(10001,*)title
	read(10001,*)wtt

	lattice_constant : do i = 1,3
	read(10001,*) (abc(i,j), j= 1,3)
	enddo lattice_constant

	allocate(ion(0:n_species),species(n_species))
	read(10001,*)(species(i), i =1,n_species)
	ion = 0
	read(10001,*)(ion(i), i =1 ,n_species)
	n_atoms = sum(ion)
	
!	allocate(x(n_atoms),y(n_atoms),z(n_atoms))

!	x = 0; y = 0; z = 0
!	read(10001,*)coordinate
!	if(coordinate.ne.'Direct')then
!	print*, 'Choose POSCAR in fraction co-ordinate only'
!	stop
!	endif
	
!	do i = 1,n_atoms
!	read(10001,*)x(i),y(i),z(i) !fraction co-ordinates
!	end do
	close(10001)

	spin_orbit : if(so.eq.'n')then

        str(1) = '# 1E-Ef'

	str(2) ='2 s_up'

	str(3) ='3 px_up'; str(4) ='4 py_up'; str(5) ='5 pz_up'
	str(6) ='6 p_up'

	str(7) ='7 dxy_up'; str(8) ='8 dyz_up'; str(9) ='9 dzx_up'
	str(10) ='10 dz2_up'; str(11) ='11 dx2_up'
	str(12) ='12 d_up'

	str(13) ='13 f-3_up'; str(14) ='14 f-2_up'; str(15) ='15 f-1_up'
	str(16) ='16 f 0_up'; str(17) ='17 f+1_up'; str(18) ='18 f+2_up'
	str(19) ='19 f+3_up'
	str(20) ='20 f_up'

	str(21) ='21 s_dn'

	str(22) ='22 px_dn'; str(23) ='23 py_dn'; str(24) ='24 pz_dn'
	str(25) ='25 p_dn'

	str(26) ='26 dxy_dn'; str(27) ='27 dyz_dn'; str(28) ='28 dzx_dn'
	str(29) ='29 dz2_dn'; str(30) ='30 dx2_dn'
	str(31) ='31 d_dn'

	str(32) ='32 f-3_dn'; str(33) ='33 f-2_dn'; str(34) ='34 f-1_dn'
	str(35) ='35 f 0_dn'; str(36) ='36 f+1_dn'; str(37) ='37 f+2_dn'
	str(38) ='38 f+3_dn'
	str(39) ='39 f_dn'
	str(40) ='40 TOTAL_up'; str(41) ='41 TOTAL_dn'
	str1(1) ='21 TOTAL'


	spin_calc: if(ispin.eq.2) then

	if (orb.eq.'s')then
	ol = 3
	orb_l = ol
	elseif(orb.eq.'p')then
	ol = 9
	orb_l = ol + 2
	elseif(orb.eq.'d')then
	ol = 19
	orb_l = ol + 4
	elseif(orb.eq.'f')then
	ol = 33
	orb_l = ol + 6
	else
	print*, 'inputs are not appropriate, please &
	chk the input file'
	stop
	end if
	
	!_____________ READING DOSCAR FILE __________!
        open(unit=10002,file='DOSCAR',status='unknown')
        read(10002,*) a1,b1,c1,d1
        read(10002,*) a2,b2,c2,d2,e2
        read(10002,*)a3
        read(10002,*)car
        read(10002,*)ttl
        read(10002,*)emax, emin, line, eferm, wt

	allocate (a(n_atoms, line, 39),b(n_atoms, line, 39),&
	 hdr(n_atoms,5), tot(5,line))

	a= 0.0; b = 0.0; hdr = 0.0

        do i= 1,line
        read(10002,*)(tot(j,line), j=1, 5)
        end do

        
	do j=1,n_atoms
        read(10002,*) (hdr(j,jj), jj = 1,5)

        do i=1,line
	read(10002,*) (a(j,i,jj), jj = 1, ol)

	b(j,i,1) = a(j,i,1); b(j,i,2) = a(j,i,2) 				! E, s_up

	b(j,i,3) = a(j,i,8) ; b(j,i,4) = a(j,i,4); b(j,i,5) = a(j,i,6) 		!px, py, pz _up
	b(j,i,6) = 0								! p_tot _up

 	b(j,i,7) = a(j,i,10); b(j,i,8) = a(j,i,12); b(j,i,9) = a(j,i,16)	!dxy, dyz, dzx _up
	b(j,i,10) = a(j,i,14); b(j,i,11) = a(j,i,18) 				!dz2, dx2 _up
	b(j,i,12) = 0								!d_tot_up	

	b(j,i,13) = a(j,i,20); b(j,i,14) = a(j,i,22); b(j,i,15) = a(j,i,24)	!f-3, f-2, f-1 _up
	b(j,i,16) = a(j,i,26); b(j,i,17) = a(j,i,28); b(j,i,18) = a(j,i,30)	!f 0, f+1, f+2 _up
	b(j,i,19) = a(j,i,32)							!f+3 _up
	b(j,i,20) = 0								!f_tot_up

	b(j,i,21) = -a(j,i,3) 							! s_dn

	b(j,i,22) = -a(j,i,9) ; b(j,i,23) = -a(j,i,5); b(j,i,24) = -a(j,i,7)	!px, py, pz _dn
	b(j,i,25) = 0								!p_tot _dn

 	b(j,i,26) = -a(j,i,11); b(j,i,27) = -a(j,i,13); b(j,i,28) = -a(j,i,17)	!dxy, dyz, dzx _dn
	b(j,i,29) = -a(j,i,15); b(j,i,30) = -a(j,i,19) 				!dz2, dx2 _dn
	b(j,i,31) = 0								!d_tot _dn

	b(j,i,32) = -a(j,i,21); b(j,i,33) = -a(j,i,23); b(j,i,34) = -a(j,i,25)	!f-3, f-2, f-1 _dn
	b(j,i,35) = -a(j,i,27); b(j,i,36) = -a(j,i,29); b(j,i,37) = -a(j,i,31)	! f 0, f+1, f+2 _dn
	b(j,i,38) = -a(j,i,33)							! f+3 _dn
	b(j,i,39) = 0								!f_tot _dn

        end do
        end do

        main1 : do j =1, n_atoms
	
	num = num + 1

        write(filename,1) num

1       format('atom_',i3.3,'.dat')

        open(file= filename,unit = num)
    
        write(num,*)'# emax =',hdr(j,1),'emin =',hdr(j,2), &
       'no. of lines =',hdr(j,3),&
        '# fermi =',hdr(j,4)!...A

	write(num,21) (str(jj),jj = 1, 41)

        do i = 1,line !do i = 2,line; if u don't want to write first line (A)

	b(j,i,6) = b(j,i,3) + b(j,i,4) + b(j,i,5)			!p_up  

	b(j,i,25) = b(j,i,22) + b(j,i,23) + b(j,i,24)			!p_dn

	b(j,i,12) = b(j,i,7) + b(j,i,8) + b(j,i,9) +&			!d_up 
	b(j,i,10) + b(j,i,11)


	b(j,i,31) = b(j,i,26) + b(j,i,27) + b(j,i,28) +&		!d_dn 
	b(j,i,29) + b(j,i,30)

	b(j,i,20) = b(j,i,13) + b(j,i,14) + b(j,i,15) +&		!f_up 
	b(j,i,16) + b(j,i,17) + b(j,i,18) + b(j,i,19) 

	b(j,i,39) = b(j,i,32) + b(j,i,33) + b(j,i,34) +&		!fdn 
	b(j,i,35) + b(j,i,36) + b(j,i,37) + b(j,i,38)

	write(num,20) b(j,i,1) - hdr(j,4), (b(j,i,jj), jj = 2, 39),tot(2,i),tot(3,i)

        end do

	close(num)

        end do main1

	atomic_tot_dos : do k = 1, n_species

        write(filename,11) species(k)

11       format('atom_',A2,'.dat')

        open(file=filename,unit= n_atoms+k)

	write(n_atoms+k,21) (str(jj),jj = 1, 41)

!	do i = int(sum(ion(0:k-1)))+1, int(sum(ion(0:k)))
	do i = 1, line
	at1 = int(sum(ion(0:k-1)))+1; at2 = int(sum(ion(0:k)))

	write(n_atoms+k,20)  b(k,i,1) - hdr(k,4), sum(b(at1:at2,i,2)), & 
	sum(b(at1:at2,i,3)),sum(b(at1:at2,i,4)),&
	sum(b(at1:at2,i,5)), sum(b(at1:at2,i,6)),&
	sum(b(at1:at2,i,7)),sum(b(at1:at2,i,8)),&
	sum(b(at1:at2,i,9)),&
	sum(b(at1:at2,i,10)), sum(b(at1:at2,i,11)), &
	sum(b(at1:at2,i,12)),&
	sum(b(at1:at2,i,13)), sum(b(at1:at2,i,14)), &
	sum(b(at1:at2,i,15)),&
	sum(b(at1:at2,i,16)), sum(b(at1:at2,i,17)), &
	sum(b(at1:at2,i,18)), sum(b(at1:at2,i,19)), &
	sum(b(at1:at2,i,20)),&
	sum(b(at1:at2,i,21)), sum(b(at1:at2,i,22)), &
	sum(b(at1:at2,i,23)), sum(b(at1:at2,i,24)), &
	sum(b(at1:at2,i,25)),&
	sum(b(at1:at2,i,26)),sum(b(at1:at2,i,27)),&
	sum(b(at1:at2,i,28)),sum(b(at1:at2,i,29)),&
	sum(b(at1:at2,i,30)),sum(b(at1:at2,i,31)),&
	sum(b(at1:at2,i,32)),&
	sum(b(at1:at2,i,33)),sum(b(at1:at2,i,34)),&
	sum(b(at1:at2,i,35)),sum(b(at1:at2,i,36)),&
	sum(b(at1:at2,i,37)),sum(b(at1:at2,i,38)),&
	sum(b(at1:at2,i,39)), &

	sum(b(at1:at2,i,2))+ & !s_up
	sum(b(at1:at2,i,6))+& !p_t_up
	sum(b(at1:at2,i,12))+& ! d_t_up
	sum(b(at1:at2,i,20)),& ! f_t_up

	sum(b(at1:at2,i,21))+ & !s_dn
	sum(b(at1:at2,i,25))+& !p_t_dn
	sum(b(at1:at2,i,31))+& ! d_t_dn
	sum(b(at1:at2,i,39)) ! f_t_dn


	end do

	close(n_atoms+k)
	end do atomic_tot_dos


	elseif(ispin.eq.1)then

	if (orb.eq.'s')then
	ol = 2
	orb_l = ol
	elseif(orb.eq.'p')then
	ol = 5
	orb_l = ol + 1
	elseif(orb.eq.'d')then
	ol = 10
	orb_l = ol + 2
	elseif(orb.eq.'f')then
	ol = 17
	orb_l = ol + 3
	else
	print*, 'inputs are not appropriate, please &
	chk the input file'
	stop
	end if

	!_____________ READING DOSCAR FILE __________!
        open(unit=10002,file='DOSCAR',status='unknown')
        read(10002,*) a1,b1,c1,d1
        read(10002,*) a2,b2,c2,d2,e2
        read(10002,*)a3
        read(10002,*)car
        read(10002,*)ttl
        read(10002,*)emax, emin, line, eferm, wt

	allocate (a(n_atoms, line, 20),b(n_atoms, line, 20),&
	 hdr(n_atoms,5), tot(3,line))

	a= 0.0; b = 0.0; hdr = 0.0

        do i= 1,line
        read(10002,*)(tot(j,line), j=1, 3)
        end do

	do j=1,n_atoms
        read(10002,*) (hdr(j,jj), jj = 1,5)

        do i=1,line
	read(10002,*) (a(j,i,jj), jj = 1, ol)

	b(j,i,1) = a(j,i,1); b(j,i,2) = a(j,i,2) 				! E, s_up

	b(j,i,3) = a(j,i,5); b(j,i,4) = a(j,i,3); b(j,i,5) = a(j,i,4) 		!px, py, pz _up
	b(j,i,6) = 0								!p_tot _up

 	b(j,i,7) = a(j,i,6); b(j,i,8) = a(j,i,7); b(j,i,9) = a(j,i,9)		!dxy, dyz, dzx _up
	b(j,i,10) = a(j,i,8); b(j,i,11) = a(j,i,10) 				!dz2, dx2 _up
	b(j,i,12) = 0								!d_tot_up

	b(j,i,13) = a(j,i,11); b(j,i,14) = a(j,i,12); b(j,i,15) = a(j,i,13)	!f-3, f-2, f-1 _up
	b(j,i,16) = a(j,i,14); b(j,i,17) = a(j,i,15); b(j,i,18) = a(j,i,16)	!f 0, f+1, f+2 _up
	b(j,i,19) = a(j,i,17)							!f+3 _up
	b(j,i,20) = 0								!f_tot_up

        end do
        end do


        main2 : do j =1, n_atoms
               num = num + 1


        write(filename,3) num

3       format('atom_',i3.3,'.dat')
!1       format('atom_',i3.3)


        open(file=filename,unit = num)
    
        write(num,*)'# emax =',hdr(j,1),'emin =',hdr(j,2), &
       'no. of lines =',hdr(j,3),&
        '# fermi =',hdr(j,4)!...A

	write(num,21) (str(jj),jj = 1, 20), str1(1)

        do i = 1,line !do i = 2,line; if u don't want to write first line (A)

	b(j,i,6) = b(j,i,3) + b(j,i,4) + b(j,i,5)			!p_up  

	b(j,i,12) = b(j,i,7) + b(j,i,8) + b(j,i,9) +&			!d_up 
	b(j,i,10) + b(j,i,11)

	b(j,i,20) = b(j,i,13) + b(j,i,14) + b(j,i,15) +&		!f_up 
	b(j,i,16) + b(j,i,17) + b(j,i,18) + b(j,i,19) 


	write(num,20) b(j,i,1) - hdr(j,4), (b(j,i,jj), jj = 2, orb_l),tot(2,i)
        end do

	close(num)

        end do main2

	atomic_tot_dos_ : do k = 1, n_species

        write(filename,12) species(k)

12       format('atom_',A2,'.dat')

        open(file=filename,unit=n_atoms+k)

	write(n_atoms+k,21) (str(jj),jj = 1, 20), str1(1)


	do i = 1, line
	at1 = int(sum(ion(0:k-1)))+1; at2 = int(sum(ion(0:k)))

	write(n_atoms+k,20)  b(k,i,1) - hdr(k,4), &
	sum(b(at1:at2,i,2)), & !s

	sum(b(at1:at2,i,3)),sum(b(at1:at2,i,4)),& !px, py
	sum(b(at1:at2,i,5)),sum(b(at1:at2,i,6)),& !pz, p_t

	sum(b(at1:at2,i,7)),sum(b(at1:at2,i,8)),& !dxy, dyz
	sum(b(at1:at2,i,9)),sum(b(at1:at2,i,10)),& !dzx, dz2
	sum(b(at1:at2,i,11)),sum(b(at1:at2,i,12)),& !dx2, d_t


	sum(b(at1:at2,i,13)),sum(b(at1:at2,i,14)),& !f-3, f-2
	sum(b(at1:at2,i,15)),sum(b(at1:at2,i,16)),& !f-1, f 0
	sum(b(at1:at2,i,17)),sum(b(at1:at2,i,18)),& !f+1 f+2
	sum(b(at1:at2,i,19)),sum(b(at1:at2,i,20)),& !f+3, f_t


	sum(b(at1:at2,i,2))+ & !s

	sum(b(at1:at2,i,6))+& !p_t

	sum(b(at1:at2,i,12))+& ! d_t

	sum(b(at1:at2,i,20)) ! f_t

	end do
	close(n_atoms+k)

	end do atomic_tot_dos_


	else
	print*,'ERROR: Select the correct ISPIN value'
	end if spin_calc
	
	elseif(so.eq.'y')then

        str(1) = '# 1E-Ef'

	str(2) ='2 s'; str(3) ='3 s_mx'; str(4) ='4 s_my'; str(5) ='5 s_mz'

	str(6) ='6 px'; str(7) ='7 px_mx'; str(8) ='8 px_my'; str(9) ='9 px_mz'
	str(10) ='10 py'; str(11) ='11 py_mx'; str(12) ='12 py_my'; str(13) ='13 py_mz'
	str(14) ='14 pz'; str(15) ='15 pz_mx'; str(16) ='16 pz_my'; str(17) ='17 pz_mz'
	str(18) ='18 p_tot'

	str(19) ='19 dxy'; str(20) ='20 dxy_mx'; str(21) ='21 dxy_my'; str(22) ='22 dxy_mz'
	str(23) ='23 dyz'; str(24) ='24 dyz_mx'; str(25) ='25 dyz_my'; str(26) ='26 dyz_mz'
	str(27) ='27 dxz'; str(28) ='28 dxz_mx'; str(29) ='29 dzx_my'; str(30) ='30 dzx_mz'
	str(31) ='31 dz2'; str(32) ='32 dz2_mx'; str(33) ='33 dz2_my'; str(34) ='34 dz2_mz'
	str(35) ='35 dx2'; str(36) ='36 dx2_mx'; str(37) ='37 dx2_my'; str(38) ='38 dx2_mz'
	str(39) ='39 d_tot'

	str(40) ='40 f-3'; str(41) ='41 f-3_mx'; str(42) ='42 f-3_my'; str(43) ='43 f-3_mz'
	str(44) ='44 f-2'; str(45) ='45 f-2_mx'; str(46) ='46 f-2_my'; str(47) ='47 f-2_mz'
	str(48) ='48 f-1'; str(49) ='49 f-1_mx'; str(50) ='50 f-1_my'; str(51) ='51 f-1_mz'
	str(52) ='52 f 0'; str(53) ='53 f 0_mx'; str(54) ='54 f 0_my'; str(55) ='55 f 0_mz'
	str(56) ='56 f+1'; str(57) ='57 f+1_mx'; str(58) ='58 f+1_my'; str(59) ='59 f+1_mz'
	str(60) ='60 f+2'; str(61) ='61 f+2_mx'; str(62) ='62 f+2_my'; str(63) ='63 f+2_mz'
	str(64) ='64 f+3'; str(65) ='65 f+3_mx'; str(66) ='66 f+3_my'; str(67) ='67 f+2_mz'
	str(68) ='68 f_tot'

	str1(1) ='69 TOTAL'

	if (orb.eq.'s')then
	ol = 5
	orb_l = ol
	elseif(orb.eq.'p')then
	ol = 17
	orb_l = ol + 1
	elseif(orb.eq.'d')then
	ol = 37
	orb_l = ol + 2
	elseif(orb.eq.'f')then
	ol = 65
	orb_l = ol + 3
	else
	print*, 'inputs are not appropriate, please &
	chk the input file'
	stop
	end if


	!_____________ READING DOSCAR FILE __________!
        open(unit=10002,file='DOSCAR',status='unknown')
        read(10002,*) a1,b1,c1,d1
        read(10002,*) a2,b2,c2,d2,e2
        read(10002,*)a3
        read(10002,*)car
        read(10002,*)ttl
        read(10002,*)emax, emin, line, eferm, wt

	allocate (a(n_atoms, line, 70),b(n_atoms, line, 70),&
	 hdr(n_atoms,5), tot(5,line))

	a= 0.0; b = 0.0; hdr = 0.0

        do i= 1,line
        read(10002,*)(tot(j,line), j=1, 3)
        end do

	do j=1,n_atoms

        read(10002,*) (hdr(j,jj), jj = 1,5)

        do i=1,line
	read(10002,*) (a(j,i,jj), jj = 1, ol)

	b(j,i,1) = a(j,i,1)
 
	b(j,i,2) = a(j,i,2); b(j,i,3) = a(j,i,3) ; b(j,i,4) = a(j,i,4)
	b(j,i,5) = a(j,i,5)				!s			

	b(j,i,6) = a(j,i,14); b(j,i,7) = a(j,i,15); b(j,i,8) = a(j,i,16)
	b(j,i,9) = a(j,i,17)				!px
	b(j,i,10) = a(j,i,6); b(j,i,11) = a(j,i,7); b(j,i,12) = a(j,i,8)
	b(j,i,13) = a(j,i,9)				!py		
	b(j,i,14) = a(j,i,10); b(j,i,15) = a(j,i,11); b(j,i,16) = a(j,i,12)
	b(j,i,17) = a(j,i,13)				!pz			
	b(j,i,18) = 0					!p_tot

	b(j,i,19) = a(j,i,18); b(j,i,20) = a(j,i,19); b(j,i,21) = a(j,i,20)
	b(j,i,22) = a(j,i,21)				!dxy
	b(j,i,23) = a(j,i,22); b(j,i,24) = a(j,i,23); b(j,i,25) = a(j,i,24)
	b(j,i,26) = a(j,i,25)				!dyz
	b(j,i,27) = a(j,i,30); b(j,i,28) = a(j,i,31); b(j,i,29) = a(j,i,32)
	b(j,i,30) = a(j,i,33)				!dzx
	b(j,i,31) = a(j,i,26); b(j,i,32) = a(j,i,27); b(j,i,33) = a(j,i,28)
	b(j,i,34) = a(j,i,29)				!dz2
	b(j,i,35) = a(j,i,34); b(j,i,36) = a(j,i,35); b(j,i,37) = a(j,i,36)
	b(j,i,38) = a(j,i,37)				!dx2 	
	b(j,i,39) = 0					!d_tot

	b(j,i,40) = a(j,i,38); b(j,i,41) = a(j,i,39); b(j,i,42) = a(j,i,40)
	b(j,i,43) = a(j,i,41)				!f-3
	b(j,i,44) = a(j,i,42); b(j,i,45) = a(j,i,43); b(j,i,46) = a(j,i,44)
	b(j,i,47) = a(j,i,45)				!f-2
	b(j,i,48) = a(j,i,46); b(j,i,49) = a(j,i,47); b(j,i,50) = a(j,i,48)
	b(j,i,51) = a(j,i,49)				!f-1
	b(j,i,52) = a(j,i,50); b(j,i,53) = a(j,i,51); b(j,i,54) = a(j,i,52)
	b(j,i,55) = a(j,i,53)				!f 0
	b(j,i,56) = a(j,i,54); b(j,i,57) = a(j,i,55); b(j,i,58) = a(j,i,56)
	b(j,i,59) = a(j,i,57)				!f+1
	b(j,i,60) = a(j,i,58); b(j,i,61) = a(j,i,59); b(j,i,62) = a(j,i,60)
	b(j,i,63) = a(j,i,61)				!f+2
	b(j,i,64) = a(j,i,62); b(j,i,65) = a(j,i,63); b(j,i,66) = a(j,i,64)
	b(j,i,67) = a(j,i,65)
	b(j,i,68) = 0					!f_tot 

        end do
        end do

        main_so : do j =1, n_atoms
	
	num = num+1

        write(filename,111) num

111     format('atom_so',i3.3,'.dat')

        open(file=filename,unit = num)
    
        write(num,*)'# emax =',hdr(j,1),'emin =',hdr(j,2), &
       'no. of lines =',hdr(j,3),&
        '# fermi =',hdr(j,4)!...A

	write(num,21) (str(jj),jj = 1, orb_l)



        do i = 1,line !do i = 2,line; if u don't want to write first line (A)

	b(j,i,18) = b(j,i,6) + b(j,i,10) + b(j,i,14)			!p_tot  


	b(j,i,39) = b(j,i,19) + b(j,i,23) + b(j,i,27) +&		!d_tot 
	b(j,i,31) + b(j,i,35)

	b(j,i,68) =  b(j,i,40) + b(j,i,44) + b(j,i,48) +&		!f_tot
	b(j,i,52) + b(j,i,56) + b(j,i,60) + b(j,i,64)

	write(num,20) b(j,i,1) - hdr(j,4), (b(j,i,jj), jj = 2, orb_l)
        end do

	close(num)

        end do main_so

	atomic_tot_dos_so : do k = 1, n_species

        write(filename,112) species(k)

112     format('atom_so',A2,'.dat')

        open(file=filename,unit= n_atoms+k)

	write(n_atoms+k,21) (str(jj),jj = 1, 68), str1(1)

	do i = 1, line
	at1 = int(sum(ion(0:k-1)))+1; at2 = int(sum(ion(0:k)))

	write(n_atoms+k,20)  b(k,i,1) - hdr(k,4), sum(b(at1:at2,i,2)), & 
	sum(b(at1:at2,i,3)),sum(b(at1:at2,i,4)),&
	sum(b(at1:at2,i,5)), sum(b(at1:at2,i,6)),&
	sum(b(at1:at2,i,7)),sum(b(at1:at2,i,8)),&
	sum(b(at1:at2,i,9)),&
	sum(b(at1:at2,i,10)), sum(b(at1:at2,i,11)), &
	sum(b(at1:at2,i,12)),&
	sum(b(at1:at2,i,13)), sum(b(at1:at2,i,14)), &
	sum(b(at1:at2,i,15)),&
	sum(b(at1:at2,i,16)), sum(b(at1:at2,i,17)), &
	sum(b(at1:at2,i,18)), sum(b(at1:at2,i,19)), &
	sum(b(at1:at2,i,20)),&
	sum(b(at1:at2,i,21)), sum(b(at1:at2,i,22)), &
	sum(b(at1:at2,i,23)), sum(b(at1:at2,i,24)), &
	sum(b(at1:at2,i,25)),&
	sum(b(at1:at2,i,26)),sum(b(at1:at2,i,27)),&
	sum(b(at1:at2,i,28)),sum(b(at1:at2,i,29)),&
	sum(b(at1:at2,i,30)),sum(b(at1:at2,i,31)),&
	sum(b(at1:at2,i,32)),&
	sum(b(at1:at2,i,33)),sum(b(at1:at2,i,34)),&
	sum(b(at1:at2,i,35)),sum(b(at1:at2,i,36)),&
	sum(b(at1:at2,i,37)),sum(b(at1:at2,i,38)),&
	sum(b(at1:at2,i,39)), &

	sum(b(at1:at2,i,40)),&
	sum(b(at1:at2,i,41)), sum(b(at1:at2,i,42)), &
	sum(b(at1:at2,i,43)),&
	sum(b(at1:at2,i,44)), sum(b(at1:at2,i,45)), &
	sum(b(at1:at2,i,46)),&
	sum(b(at1:at2,i,47)), sum(b(at1:at2,i,48)), &
	sum(b(at1:at2,i,49)), sum(b(at1:at2,i,50)), &
	sum(b(at1:at2,i,51)),&
	sum(b(at1:at2,i,52)), sum(b(at1:at2,i,53)), &
	sum(b(at1:at2,i,54)), sum(b(at1:at2,i,55)), &
	sum(b(at1:at2,i,56)),&
	sum(b(at1:at2,i,57)),sum(b(at1:at2,i,58)),&
	sum(b(at1:at2,i,59)),sum(b(at1:at2,i,60)),&
	sum(b(at1:at2,i,61)),sum(b(at1:at2,i,62)),&
	sum(b(at1:at2,i,63)),&
	sum(b(at1:at2,i,64)),sum(b(at1:at2,i,65)),&
	sum(b(at1:at2,i,66)),sum(b(at1:at2,i,67)),&
	sum(b(at1:at2,i,68)),&

	sum(b(at1:at2,i,2))+ & !s_tot
	sum(b(at1:at2,i,18))+& !p_tot
	sum(b(at1:at2,i,39))+& ! d_tot
	sum(b(at1:at2,i,68)) ! f_tot


	end do

	close(n_atoms+k)

	end do atomic_tot_dos_so


	else

	print*, "Please specify correct 'SO' information in input file"
	stop

	end if spin_orbit

	call system('mkdir plot_files')
	call system ('bash reduce.sh')
        call system('mv atom* plot_files')
20	format(70f12.7)
21      format(70A12)

        stop
        end
