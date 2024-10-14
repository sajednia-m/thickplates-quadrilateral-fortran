!=====================================================================
!=====================================================================
!======================= FINITE ELEMENT METHOD =======================
!====== THICK PLATES USING QUADRILATERAL ELEMENTS WITH 8 NODES =======
!========================= MOHAMMAD SAJEDNIA =========================
!=====================================================================
!=====================================================================
program thickplate_q8
  implicit none
  integer :: dimen, nne, nodof, eldof, nnl, nbf, nxe, nye, ngpb, ngps
  integer :: i, nel, nnd, nabc
  integer, dimension(:,:), allocatable :: elements, u_cnt
  real(8), dimension(:,:), allocatable :: bc, add_bc, nodes, nloads
  real(8), dimension(:,:), allocatable :: bodyforces, kk, f, u
  real(8), dimension(:,:), allocatable :: rf, ef, nforce
  real(8), dimension(3,3) :: deeb
  real(8), dimension(2,2) :: dees
  real(8) :: E, nu, thick, length, width
  real(8) :: dhx, dhy
  character(len=20) :: filename

!=====================================================================
!========================== READING INPUTS ===========================
!=====================================================================

  filename = 'inputs.txt'
  open(unit=10, file=filename, status='old', action='read')

  do i=1,4
    read(10,*)
  end do
  read(10,*) dimen, thick

  read(10,*)
  read(10,*)
  read(10,*) nne, nodof
  eldof = nne * nodof ! 24=8*3

  read(10,*)
  read(10,*)
  read(10,*) length, width
  read(10,*)
  read(10,*)
  read(10,*) nxe, nye
  dhx = length / nxe
  dhy = width / nye

  do i=1,4
    read(10,*)
  end do
  read(10,*) e, nu
  do i=1,7
    read(10,*)
  end do

  allocate(bc(4,3))
  read(10,*) bc(1,:)
  read(10,*)
  read(10,*)
  read(10,*) bc(2,:)
  read(10,*)
  read(10,*)
  read(10,*) bc(3,:)
  read(10,*)
  read(10,*)
  read(10,*) bc(4,:)
  read(10,*)
  read(10,*) nabc
  read(10,*)
  if (nabc /= 0) then
    allocate(add_bc(nabc,4))
    do i=1,nabc
        read(10,*) add_bc(i,:)
    end do
    else
        read(10,*)
  end if
  do i=1,5
    read(10,*)
  end do

  read(10,*) nnl
  read(10,*)
  if (nnl /= 0) then
    allocate(nloads(nnl,4))
    do i=1,nnl
        read(10,*) nloads(i,:)
    end do
    else
        read(10,*)
  end if
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*) nbf
  read(10,*)
  if (nbf /= 0) then
    allocate(bodyforces(nbf,2))
    do i=1,nbf
        read(10,*) bodyforces(i,:)
    end do
    else
        read(10,*)
  end if

  do i=1,4
    read(10,*)
  end do
  read(10,*) ngpb, ngps

  nel = nxe*nye
  nnd = (1+2*nxe)*(1+2*nye)-(nxe*nye)
  allocate(nodes(nnd,2))
  allocate(elements(nel,8))
  call q8_mesh(length,width,nxe,nye,dhx,dhy,nnd,nel,nodes,elements)

  call formdeeb(e,nu,thick,deeb)
  call formdees(e,nu,thick,dees)

  allocate(f(3*nnd,1))
  call forces(ngpb,nnd,nodes,nel,elements,nnl,nloads,nbf,bodyforces,f)

  allocate(kk(nnd*3,nnd*3))
  call stiffness(ngpb,ngps,thick,nodes,elements,nel,nnd,deeb,dees,kk)

  allocate(u(3*nnd,1))
  allocate(u_cnt(3*nnd,1))
  call displacements(length, width, nnd, nabc, bc, add_bc, nodes, kk, f, u, u_cnt)

  allocate(rf(3*nnd,1))
  call reactions(nnd, kk, u_cnt, u, rf)

  allocate(ef(nel,5))
  call element_forces(nodes,elements,nel,nnd,deeb,dees,u,ef)

  allocate(nforce(nnd,5))
  call nodes_forces(elements,nel,nnd,ef,nforce)

!=====================================================================
!========================== WRITING RESULTS ==========================
!=====================================================================

  filename = 'outputs.txt'
  open(unit=16, file=filename, status='replace', action='write')

  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '==================================== DISPLACEMENTS ===================================='
  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '        Node     X_Coor      Y_Coor      W        ThetaX      ThetaY'
  write(16, '(A)') '---------------------------------------------------------------------------------------'
  do i = 1, nnd
      write(16, '(I10, 2F12.2, 3F12.6)') i, nodes(i,1), nodes(i,2), u(3*i-2,1), u(3*i-1,1), u(3*i,1)
  end do
  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '=================================== REACTION FORCES ==================================='
  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '        Node        W             ThetaX        ThetaY'
  write(16, '(A)') '---------------------------------------------------------------------------------------'
  do i = 1, nnd
      write(16, '(I10, 3F15.4)') i, rf(3*i-2,1), rf(3*i-1,1), rf(3*i,1)
  end do
  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '================================= MOMENTS AND SHEARS =================================='
  write(16, '(A)') '======================================================================================='
  write(16, '(A)') '        Node       Mx             My             Mxy            Qx            Qy'
  write(16, '(A)') '---------------------------------------------------------------------------------------'
  do i = 1, nnd
      write(16, '(I10, 5F15.4)') i, nforce(i,1), nforce(i,2), nforce(i,3), nforce(i,4), nforce(i,5)
  end do
  write(16, '(A)') '======================================================================================='
  close(16)
  print*, '==========================================================='
  print*, '==================== Check outputs.txt ===================='
  print*, '==========================================================='

  deallocate(nodes, elements, bc, u, kk, u_cnt, rf, ef, nforce)

  if (nabc /= 0) then
    deallocate(add_bc)
  end if
  if (nnl /= 0) then
    deallocate(nloads)
  end if
  if (nbf /= 0) then
    deallocate(bodyforces)
  end if

  stop

contains

!=====================================================================
!===================== ASSEMBLING FORCES VECTOR ======================
!=====================================================================

subroutine forces(ngpb,nnd,nodes,nel,elements,nnl,nloads,nbf,bodyforces,f)
    implicit none
    integer, intent(in) :: ngpb, nnd, nel, nnl, nbf
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(:,:), intent(in) :: nodes, nloads, bodyforces
    real(8), dimension(3*nnd,1), intent(out) :: f
    integer :: i, j, ig, jg, m, nload_dof, dof, bf_dof
    real(8) :: nload_x, nload_y, nload_mag, x, y, wi, wj, bf_mag, detj
    real(8), dimension(ngpb,2) :: sampb
    real(8), dimension(8,2) :: eng
    real(8), dimension(8,1) :: fun
    real(8), dimension(2,8) :: der
    real(8), dimension(2,2) :: dj, dji
    integer, dimension(8,1) :: n
    f = 0.0

    if (nnl/=0) then
        do i=1,nnl
            nload_x = nloads(i,1)
            nload_y = nloads(i,2)
            nload_mag = nloads(i,3)
            nload_dof = nloads(i,4)
            do j=1,nnd
                x = nodes(j,1)
                y = nodes(j,2)
                if (abs(nload_x - x)<=0.00001 .AND. abs(nload_y - y)<=0.00001) then
                    f(3*(j-1)+nload_dof,1) = f(3*(j-1)+nload_dof,1) + nload_mag
                end if
            end do
        end do
    end if

    call gauss_points(ngpb,sampb)

    if (nbf/=0) then
        do i=1,nbf
            bf_mag = bodyforces(i,1)
            bf_dof = bodyforces(i,2)
            do j=1,nel
                call element(j,nodes,elements,eng)
                n(1,1) = elements(j,1)
                n(2,1) = elements(j,2)
                n(3,1) = elements(j,3)
                n(4,1) = elements(j,4)
                n(5,1) = elements(j,5)
                n(6,1) = elements(j,6)
                n(7,1) = elements(j,7)
                n(8,1) = elements(j,8)
                do ig=1,ngpb
                    wi=sampb(ig,2)
                    do jg=1,ngpb
                        wj=sampb(jg,2)
                        call shape_q8(sampb,ig,jg,fun,der)
                        call sub_matmul(der, eng, dj, 2, 8, 2)
                        call jac_inv(dj, dji, detj)
                        do m=1,8
                            f(3*(n(m,1)-1)+bf_dof,1) = f(3*(n(m,1)-1)+bf_dof,1) + detj * wi * wj * thick * bf_mag * fun(m,1)
                        end do
                    end do
                end do
            end do
        end do
    end if

end subroutine forces

!=====================================================================
!========== MESH USING QUADRILATERAL ELEMENTS WITH 8 NODES ===========
!=====================================================================

subroutine q8_mesh(length,width,nxe,nye,dhx,dhy,nnd,nel,nodes,elements)
    implicit none
    integer, intent(in) :: nxe, nye
    real(8), intent(in) :: length, width, dhx, dhy
    integer, intent(in) :: nnd, nel
    integer, dimension(nel,8), intent(out) :: elements
    real(8), dimension(nnd,2), intent(out) :: nodes
    integer :: i,j,k,n1,n2,n3,n4,n5,n6,n7,n8
    real(8) :: x_origin, y_origin

    nodes = 0.0
    elements = 0.0
    x_origin = 0.0
    y_origin = 0.0
    k = 0

    do i=1,nxe
        do j=1,nye
            k=k+1
            n1 = (i-1)*(3*nye+2)+2*j - 1
            n2 = i*(3*nye+2)+j - nye - 1
            n3 = i*(3*nye+2)+2*j-1
            n4 = n3 + 1
            n5 = n3 + 2
            n6 = n2 + 1
            n7 = n1 + 2
            n8 = n1 + 1
            nodes(n1,1) = (i-1)*dhx - x_origin
            nodes(n1,2) = (j-1)*dhy - y_origin
            nodes(n3,1) = i*dhx - x_origin
            nodes(n3,2) = (j-1)*dhy - y_origin
            nodes(n2,1) = (nodes(n1,1)+nodes(n3,1))/2
            nodes(n2,2) = (nodes(n1,2)+nodes(n3,2))/2
            nodes(n5,1) = i*dhx- x_origin
            nodes(n5,2) = j*dhy - y_origin
            nodes(n4,1) = (nodes(n3,1)+ nodes(n5,1))/2
            nodes(n4,2) = (nodes(n3,2)+ nodes(n5,2))/2
            nodes(n7,1) = (i-1)*dhx - x_origin
            nodes(n7,2) = j*dhy - y_origin
            nodes(n6,1) = (nodes(n5,1)+ nodes(n7,1))/2
            nodes(n6,2) = (nodes(n5,2)+ nodes(n7,2))/2
            nodes(n8,1) = (nodes(n1,1)+ nodes(n7,1))/2
            nodes(n8,2) = (nodes(n1,2)+ nodes(n7,2))/2
            elements(k,1) = n1
            elements(k,2) = n2
            elements(k,3) = n3
            elements(k,4) = n4
            elements(k,5) = n5
            elements(k,6) = n6
            elements(k,7) = n7
            elements(k,8) = n8
        end do
    end do

end subroutine q8_mesh

!=====================================================================
!========================== FORM DB MATRIX ===========================
!=====================================================================

subroutine formdeeb(e,nu,thick,deeb)
    implicit none
    real(8), intent(in) :: e, nu, thick
    real(8), dimension(3,3), intent(out) :: deeb
    real(8) :: dr

    dr = e * thick**3 / (12*(1-nu*nu))
    deeb(1,1) = dr
    deeb(1,2) = dr*nu
    deeb(1,3) = 0.0
    deeb(2,1) = dr*nu
    deeb(2,2) = dr
    deeb(2,3) = 0.0
    deeb(3,1) = 0.0
    deeb(3,2) = 0.0
    deeb(3,3) = dr*(1-nu)/2

end subroutine formdeeb

!=====================================================================
!========================== FORM DS MATRIX ===========================
!=====================================================================

subroutine formdees(e,nu,thick,dees)
    implicit none
    real(8), intent(in) :: e, nu, thick
    real(8), dimension(2,2), intent(out) :: dees
    real(8) :: g

    g = e / (2*(1+nu))
    dees(1,1) = g*thick
    dees(1,2) = 0.0
    dees(2,1) = 0.0
    dees(2,2) = g*thick

end subroutine formdees

!=====================================================================
!============= APPLYING BCS AND CALCULATE DISPLACEMENTS ==============
!=====================================================================

subroutine displacements(length, width, nnd, nabc, bc, add_bc, nodes, k, f, u, u_cnt)
    implicit none
    integer, intent(in) :: nnd, nabc
    real(8), intent(in) :: length, width
    real(8), dimension(:,:), intent(in) :: nodes, bc, add_bc, k
    real(8), dimension(3*nnd,1), intent(out) :: u
    integer :: i, j, bc_cnt, dof, adof, cnti, cntj, cntii, cntf
    integer :: cond1, cond2, cond3
    real(8) :: left1, left2, left3, bottom1, bottom2, bottom3, right1, right2, right3, up1, up2, up3
    real(8) :: x, y, x_cor, y_cor, valu
    integer, dimension(3*nnd,1) :: u_cnt
    real(8), dimension(3*nnd,1) :: f
    real(8), allocatable, dimension(:) :: f_f, uf
    real(8), allocatable, dimension(:,:) :: k_f

    u = 0.0
    u_cnt = 0.0
    bc_cnt = 0

    left1 = bc(1,1)
    left2 = bc(1,2)
    left3 = bc(1,3)
    bottom1 = bc(2,1)
    bottom2 = bc(2,2)
    bottom3 = bc(2,3)
    right1 = bc(3,1)
    right2 = bc(3,2)
    right3 = bc(3,3)
    up1 = bc(4,1)
    up2 = bc(4,2)
    up3 = bc(4,3)
    do i=1,nnd
        cond1 = 0
        cond2 = 0
        cond3 = 0
        if (nodes(i,1)==0) then
            if(left1/=1000)then
                u(3*(i-1)+1,1) = left1
                f(:,1) = f(:,1) - k(:,3*(i-1)+1)*left1
                bc_cnt = bc_cnt + 1
                u_cnt(3*(i-1)+1,1) = bc_cnt
                cond1 = 1
            end if
            if(left2/=1000)then
                u(3*(i-1)+2,1) = left2
                f(:,1) = f(:,1) - k(:,3*(i-1)+2)*left2
                bc_cnt = bc_cnt + 1
                u_cnt(3*(i-1)+2,1) = bc_cnt
                cond2 = 1
            end if
            if(left3/=1000)then
                u(3*(i-1)+3,1) = left3
                f(:,1) = f(:,1) - k(:,3*(i-1)+3)*left3
                bc_cnt = bc_cnt + 1
                u_cnt(3*(i-1)+3,1) = bc_cnt
                cond3 = 1
            end if
        end if
        if (nodes(i,2)==0) then
            if(bottom1/=1000)then
                if (cond1==0) then
                    u(3*(i-1)+1,1) = bottom1
                    f(:,1) = f(:,1) - k(:,3*(i-1)+1)*bottom1
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+1,1) = bc_cnt
                    cond1 = 1
                end if
            end if
            if(bottom2/=1000)then
                if (cond2==0) then
                    u(3*(i-1)+2,1) = bottom2
                    f(:,1) = f(:,1) - k(:,3*(i-1)+2)*bottom2
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+2,1) = bc_cnt
                    cond2 = 1
                end if
            end if
            if(bottom3/=1000)then
                if (cond3==0) then
                    u(3*(i-1)+3,1) = bottom3
                    f(:,1) = f(:,1) - k(:,3*(i-1)+3)*bottom3
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+3,1) = bc_cnt
                    cond3 = 1
                end if
            end if
        end if
        if (nodes(i,1)==length) then
            if(right1/=1000)then
                if (cond1==0) then
                    u(3*(i-1)+1,1) = right1
                    f(:,1) = f(:,1) - k(:,3*(i-1)+1)*right1
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+1,1) = bc_cnt
                    cond1 = 1
                end if
            end if
            if(right2/=1000)then
                if (cond2==0) then
                    u(3*(i-1)+2,1) = right2
                    f(:,1) = f(:,1) - k(:,3*(i-1)+2)*right2
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+2,1) = bc_cnt
                    cond2 = 1
                end if
            end if
            if(right3/=1000)then
                if (cond3==0) then
                    u(3*(i-1)+3,1) = right3
                    f(:,1) = f(:,1) - k(:,3*(i-1)+3)*right3
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+3,1) = bc_cnt
                    cond3 = 1
                end if
            end if
        end if
        if (nodes(i,2)==width) then
            if(up1/=1000)then
                if (cond1==0) then
                    u(3*(i-1)+1,1) = up1
                    f(:,1) = f(:,1) - k(:,3*(i-1)+1)*up1
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+1,1) = bc_cnt
                    cond1 = 1
                end if
            end if
            if(up2/=1000)then
                if (cond2==0) then
                    u(3*(i-1)+2,1) = up2
                    f(:,1) = f(:,1) - k(:,3*(i-1)+2)*up2
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+2,1) = bc_cnt
                    cond2 = 1
                end if
            end if
            if(up3/=1000)then
                if (cond3==0) then
                    u(3*(i-1)+3,1) = up3
                    f(:,1) = f(:,1) - k(:,3*(i-1)+3)*up3
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(i-1)+3,1) = bc_cnt
                    cond3 = 1
                end if
            end if
        end if
    end do

    if (nabc/=0) then
        do i=1,nabc
            x_cor = add_bc(i,1)
            y_cor = add_bc(i,2)
            valu = add_bc(i,3)
            dof = add_bc(i,4)
            do j=1,nnd
                x = nodes(j,1)
                y = nodes(j,2)
                if (abs(x_cor-x)<=0.00001 .AND. abs(y_cor - y)<=0.00001) then
                    u(3*(j-1)+dof,1) = valu
                    f(:,1) = f(:,1) - k(:,3*(j-1)+dof)*valu
                    bc_cnt = bc_cnt + 1
                    u_cnt(3*(j-1)+dof,1) = bc_cnt
                end if
            end do
        end do
    end if

    adof = 3*nnd - bc_cnt ! Active Degrees Of Freedom
    allocate(k_f(adof,adof))
    k_f = 0.0
    cnti = 1
    do i=1,3*nnd
        if (u_cnt(i,1)==0) then
            cntj = 1
            do j=1,3*nnd
                if (u_cnt(j,1)==0) then
                    k_f(cnti,cntj) = k(i,j)
                    cntj=cntj+1
                end if
            end do
            cnti=cnti+1
        end if
    end do

    allocate(f_f(adof))
    f_f = 0.0
    cntii = 1
    do i=1,3*nnd
        if (u_cnt(i,1)==0) then
            f_f(cntii) = f(i,1)
            cntii = cntii + 1
        end if
    end do

    allocate(uf(adof))
    call gauss_solver(k_f, f_f, uf)

    cntf = 1
    do i=1,3*nnd
        if (u_cnt(i,1)==0) then
            u(i,1) = uf(cntf)
            cntf=cntf+1
        end if
    end do

    deallocate(k_f, f_f, uf)

end subroutine displacements

!=====================================================================
!========================== REACTION FORCE ===========================
!=====================================================================

subroutine reactions(nnd, kk, u_cnt, u, rf)
    implicit none
    integer, intent(in) :: nnd
    integer, dimension(:,:), intent(in) :: u_cnt
    real(8), dimension(:,:), intent(in) :: kk, u
    real(8), dimension(3*nnd,1), intent(out) :: rf
    integer :: i, j
    rf = 0.0

    do i=1,3*nnd
        if (u_cnt(i,1)/=0) then
            call sub_matmul(kk(i,:), u, rf(i,1), 1, 3*nnd, 1)
        end if
    end do

end subroutine reactions

!=====================================================================
!========================= STIFFNESS MATRIX ==========================
!=====================================================================

subroutine stiffness(ngpb,ngps,thick,nodes,elements,nel,nnd,deeb,dees,kk)
    implicit none
    integer, intent(in) :: ngpb, ngps, nel, nnd
    real(8), intent(in) :: thick
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(:,:), intent(in) :: nodes, deeb, dees
    real(8), dimension(nnd*3,nnd*3), intent(out) :: kk
    integer :: i, j, ig, jg
    real(8) :: wi, wj, detj
    real(8), dimension(ngpb,2) :: sampb
    real(8), dimension(ngps,2) :: samps
    real(8), dimension(8,2) :: eng
    real(8), dimension(24,24) :: keb, kes, k1b, k1s
    real(8), dimension(8,1) :: fun
    real(8), dimension(2,8) :: der, derxy
    real(8), dimension(2,2) :: dj, dji
    real(8), dimension(3,24) :: beeb
    real(8), dimension(2,24) :: bees
    real(8), dimension(24,3) :: beeb_tr, bbtr_db
    real(8), dimension(24,2) :: bees_tr, bstr_ds

    call gauss_points(ngpb,sampb)
    call gauss_points(ngps,samps)

    kk = 0.0
    do i=1,nel
        ! Element Nodes Geometry
        call element(i,nodes,elements,eng)
        keb = 0.0
        kes = 0.0

        do ig=1,ngpb
            wi=sampb(ig,2)
            do jg=1,ngpb
                wj=sampb(jg,2)

                call shape_q8(sampb,ig,jg,fun,der)

                call sub_matmul(der, eng, dj, 2, 8, 2)

                call jac_inv(dj, dji, detj)

                call sub_matmul(dji, der, derxy, 2, 2, 8)

                call formbeeb(derxy,beeb)

                call sub_transpose(beeb, beeb_tr)
                call sub_matmul(beeb_tr, deeb, bbtr_db, 24, 3, 3)
                call sub_matmul(bbtr_db, beeb, k1b, 24, 3, 24)

                keb = keb + detj * wi * wj * k1b
            end do
        end do

        call assemble_kk(nnd, keb, i, elements, kk)

        do ig=1,ngps
            wi=samps(ig,2)
            do jg=1,ngps
                wj=samps(jg,2)

                call shape_q8(samps,ig,jg,fun,der)

                call sub_matmul(der, eng, dj, 2, 8, 2)

                call jac_inv(dj, dji, detj)

                call sub_matmul(dji, der, derxy, 2, 2, 8)

                call formbees(fun,derxy,bees)

                call sub_transpose(bees, bees_tr)
                call sub_matmul(bees_tr, dees, bstr_ds, 24, 2, 2)
                call sub_matmul(bstr_ds, bees, k1s, 24, 2, 24)

                kes = kes + detj * wi * wj * k1s * 5/6
            end do
        end do

        call assemble_kk(nnd, kes, i, elements, kk)

    end do

end subroutine stiffness

!=====================================================================
!======================== ASSEMBLE KK MATRIX =========================
!=====================================================================

subroutine assemble_kk(nnd, ke, i, elements, kk)
    implicit none
    integer, intent(in) :: nnd, i
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(24,24), intent(in) :: ke
    real(8), dimension(nnd*3,nnd*3) :: kk
    integer :: j, l, k, p
    integer, dimension(8,1) :: n

    n(1,1) = elements(i,1)
    n(2,1) = elements(i,2)
    n(3,1) = elements(i,3)
    n(4,1) = elements(i,4)
    n(5,1) = elements(i,5)
    n(6,1) = elements(i,6)
    n(7,1) = elements(i,7)
    n(8,1) = elements(i,8)

    do j=1,8
        do k=1,8
            do l=1,3
                do p=1,3
                kk(3*(n(j,1)-1)+l , 3*(n(k,1)-1)+p) = kk(3*(n(j,1)-1)+l , 3*(n(k,1)-1)+p) + ke(3*(j-1)+l , 3*(k-1)+p)
                end do
            end do
        end do
    end do

end subroutine assemble_kk

!=====================================================================
!===================== ELEMENT FORCES AT CENTER ======================
!=====================================================================

subroutine element_forces(nodes,elements,nel,nnd,deeb,dees,u,ef)
    implicit none
    integer, intent(in) :: nel, nnd
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(:,:), intent(in) :: nodes, deeb, dees, u
    real(8), dimension(nel,5), intent(out) :: ef
    integer :: i, j, l, ig, jg, ngp, node
    real(8) :: wi, wj, detj
    real(8), dimension(1,2) :: samp
    real(8), dimension(8,2) :: eng
    real(8), dimension(24,24) :: keb, kes, k1b, k1s
    real(8), dimension(8,1) :: fun
    real(8), dimension(2,8) :: der, derxy
    real(8), dimension(2,2) :: dj, dji
    real(8), dimension(3,24) :: beeb
    real(8), dimension(2,24) :: bees
    real(8), dimension(24,3) :: beeb_tr, bbtr_db
    real(8), dimension(24,2) :: bees_tr, bstr_ds
    real(8), dimension(24,1) :: eld
    real(8), dimension(3,1) :: cur_b, moment
    real(8), dimension(2,1) :: cur_s, shear
    ngp = 1
    call gauss_points(ngp,samp)

    do i=1,nel
        ! Element Nodes Geometry
        call element(i,nodes,elements,eng)
        eld = 0.0

        do j=1,8
            do l=1,3
                node = elements(i,j)
                eld(3*(j-1)+l,1)=u(3*(node-1)+l,1)
            end do
        end do

        do ig=1,ngp
            wi=samp(ig,2)
            do jg=1,ngp
                wj=samp(jg,2)
                call shape_q8(samp,ig,jg,fun,der)
                call sub_matmul(der, eng, dj, 2, 8, 2)
                call jac_inv(dj, dji, detj)
                call sub_matmul(dji, der, derxy, 2, 2, 8)
                call formbeeb(derxy,beeb)
                call formbees(fun,derxy,bees)

                call sub_matmul(beeb,eld,cur_b,3,24,1) ! Bending Curvature
                call sub_matmul(deeb,cur_b,moment,3,3,1) ! Moments
                call sub_matmul(bees,eld,cur_s,2,24,1) ! Shear Curvatures
                call sub_matmul(dees,cur_s,shear,2,2,1) ! Shear Forces

                ef(i,1) = ef(i,1) + moment(1,1)
                ef(i,2) = ef(i,2) + moment(2,1)
                ef(i,3) = ef(i,3) + moment(3,1)
                ef(i,4) = ef(i,4) + shear(1,1)
                ef(i,5) = ef(i,5) + shear(2,1)
            end do
        end do
    end do

end subroutine element_forces

!=====================================================================
!=========================== NODES FORCES ============================
!=====================================================================

subroutine nodes_forces(elements,nel,nnd,ef,nforce)
    implicit none
    integer, intent(in) :: nel, nnd
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(:,:), intent(in) :: ef
    real(8), dimension(nnd,5), intent(out) :: nforce
    integer :: i, j, k, ne
    real(8) :: mx, my, mxy, qx, qy
    nforce = 0.0

    do k=1,nnd
        mx = 0.0
        my = 0.0
        mxy = 0.0
        qx = 0.0
        qy = 0.0
        ne = 0
        do i=1,nel
            do j=1,8
                if (elements(i,j)==k) then
                    ne = ne + 1
                    mx = mx + ef(i,1)
                    my = my + ef(i,2)
                    mxy = mxy + ef(i,3)
                    qx = qx + ef(i,4)
                    qy = qy + ef(i,5)
                end if
            end do
        end do
        nforce(k,1) = mx/ne
        nforce(k,2) = my/ne
        nforce(k,3) = mxy/ne
        nforce(k,4) = qx/ne
        nforce(k,5) = qy/ne
    end do

end subroutine nodes_forces

!=====================================================================
!========================== FORM BB MATRIX ===========================
!=====================================================================

subroutine formbeeb(derxy,beeb)
    implicit none
    real(8), dimension(2,8), intent(in) :: derxy
    real(8), dimension(3,24), intent(out) :: beeb
    integer :: m, k, j
    real(8) :: x, y

    beeb = 0.0
    do m=1,8
        k = 3*m
        j = k-1
        x = derxy(1,m)
        beeb(1,j) = x
        beeb(3,k) = x
        y = derxy(2,m)
        beeb(2,k) = y
        beeb(3,j) = y
    end do

end subroutine formbeeb

!=====================================================================
!========================== FORM BS MATRIX ===========================
!=====================================================================

subroutine formbees(fun,derxy,bees)
    implicit none
    real(8), dimension(8,1), intent(in) :: fun
    real(8), dimension(2,8), intent(in) :: derxy
    real(8), dimension(2,24), intent(out) :: bees
    integer :: m, k, j, i
    real(8) :: x, y

    bees = 0.0
    do m=1,8
        k = 3*m
        j = k-1
        i = k-2
        x = derxy(1,m)
        y = derxy(2,m)
        bees(2,i) = -x
        bees(1,i) = -y
        bees(1,k) = fun(m,1)
        bees(2,j) = fun(m,1)
    end do

end subroutine formbees

!=====================================================================
!============ DETERMINATE AND INVERSE OF JACOBIAN MATRIX =============
!=====================================================================

subroutine jac_inv(dj, dji, detj)
    implicit none
    real(8), dimension(2,2), intent(in) :: dj
    real(8), dimension(2,2), intent(out) :: dji
    real(8), intent(out) :: detj

    detj = dj(1,1)*dj(2,2) - dj(1,2)*dj(2,1)
    dji(1,1) = dj(2,2)/detj
    dji(2,2) = dj(1,1)/detj
    dji(1,2) = -dj(1,2)/detj
    dji(2,1) = -dj(2,1)/detj

end subroutine jac_inv

!=====================================================================
!================ Q8 SHAPE FUNCTIONS AND DERIVATIONS =================
!=====================================================================

subroutine shape_q8(samp,ig,jg,fun,der)
    implicit none
    integer, intent(in) :: ig, jg
    real(8), dimension(:,:), intent(in) :: samp
    real(8), dimension(8,1), intent(out) :: fun
    real(8), dimension(2,8), intent(out) :: der
    real(8) :: xi, xim, xip, eta, etam, etap

    xi = samp(ig,1)
    xim = 1.0-xi
    xip = 1.0+xi
    eta = samp(jg,1)
    etam = 1.0-eta
    etap = 1.0+eta

    fun(1,1) = -0.25*xim*etam*(1.0+xi+eta)
    fun(2,1) = 0.5*etam*(1.0-xi*xi)
    fun(3,1) = -0.25*xip*etam*(1.0-xi+eta)
    fun(4,1) = 0.5*xip*(1.0-eta*eta)
    fun(5,1) = -0.25*xip*etap*(1.0-xi-eta)
    fun(6,1) = 0.5*etap*(1.0-xi*xi)
    fun(7,1) = -0.25*xim*etap*(1.0+xi-eta)
    fun(8,1) = 0.5*xim*(1.0-eta*eta)

    der(1,1)=0.25*etam*(2.0*xi + eta)
    der(1,2)=-1.0*etam*xi
    der(1,3)=0.25*etam*(2.0*xi-eta)
    der(1,4)=0.5*(1-eta*eta)
    der(1,5)=0.25*etap*(2.0*xi+eta)
    der(1,6)=-1.0*etap*xi
    der(1,7)=0.25*etap*(2.0*xi-eta)
    der(1,8)=-0.5*(1.0-eta*eta)
    der(2,1)=0.25*xim*(2.0*eta+xi)
    der(2,2)=-0.5*(1. - xi*xi)
    der(2,3)=-0.25*xip*(xi-2.0*eta)
    der(2,4)=-1.0*xip*eta
    der(2,5)=0.25*xip*(xi+2.0*eta)
    der(2,6)=0.5*(1.0-xi*xi)
    der(2,7)=-0.25*xim*(xi-2.0*eta)
    der(2,8)=-1.0*xim*eta

end subroutine shape_q8

!=====================================================================
!========================== ELEMENTS NODES ===========================
!=====================================================================

subroutine element(i,nodes,elements,eng)
    implicit none
    integer, intent(in) :: i
    integer, dimension(:,:), intent(in) :: elements
    real(8), dimension(:,:), intent(in) :: nodes
    real(8), dimension(8,2), intent(out) :: eng
    integer :: j, k, l, current_node

    eng = 0.0
    do k=1,8
        do j=1,2
            current_node = elements(i,k)
            eng(k,j) = nodes(current_node,j)
        end do
    end do

end subroutine element

!=====================================================================
!=========================== GAUSS POINTS ============================
!=====================================================================

subroutine gauss_points(ngp,samp)
    implicit none
    integer, intent(in) :: ngp
    real(8), dimension(ngp,2), intent(out) :: samp

    samp = 0.0
    if (ngp==1) then
        samp(1,1)=0.0
        samp(1,2)=2.0
    end if

    if (ngp==2) then
        samp(1,1)=-1.0/sqrt(3.0)
        samp(1,2)=1.0
        samp(2,1)=1.0/sqrt(3.0)
        samp(2,2)=1.0
    end if

    if (ngp==3) then
        samp(1,1)=-0.2*sqrt(15.0)
        samp(1,2)=5.0/9.0
        samp(2,1)=0.0
        samp(2,2)=8.0/9.0
        samp(3,1)=0.2*sqrt(15.0)
        samp(3,2)=5.0/9.0
    end if

    if (ngp==4) then
        samp(1,1)=-0.861136311594053
        samp(1,2)=0.347854845137454
        samp(2,1)=-0.339981043584856
        samp(2,2)=0.652145154862546
        samp(3,1)=0.339981043584856
        samp(3,2)=0.652145154862546
        samp(4,1)=0.861136311594053
        samp(4,2)=0.347854845137454
    end if

end subroutine gauss_points

!=====================================================================
!==================== GAUSSIAN ELIMINATION METHOD ====================
!=====================================================================

subroutine gauss_solver(k,f,u)
  implicit none
  integer :: n
  real(8), dimension(:), intent(inout) :: u, f
  real(8), dimension(:,:), intent(inout) :: k
  real :: temp
  integer :: i, j, p, q
  real :: factor
  n = size(u)

  do i = 1, n-1
    p = i
    do j = i+1, n
      if (abs(k(j,i)) > abs(k(p,i))) then
        p = j
      end if
    end do

    if (p /= i) then
      do q = 1, n
        temp = k(i,q)
        k(i,q) = k(p,q)
        k(p,q) = temp
      end do
      temp = f(i)
      f(i) = f(p)
      f(p) = temp
    end if

    do j = i+1, n
      factor = k(j,i) / k(i,i)
      k(j,i) = 0.0
      do q = i+1, n
        k(j,q) = k(j,q) - factor*k(i,q)
      end do
      f(j) = f(j) - factor*f(i)
    end do
  end do

  u(n) = f(n) / k(n,n)
  do i = n-1, 1, -1
    u(i) = f(i)
    do j = i+1, n
      u(i) = u(i) - k(i,j)*u(j)
    end do
    u(i) = u(i) / k(i,i)
  end do

end subroutine gauss_solver

!=====================================================================
!===================== MATRICE'S MULTIPLICATION ======================
!=====================================================================

subroutine sub_matmul(A, B, C, m, n, p)
  implicit none
  integer :: m, n, p
  real(8) :: A(m,n), B(n,p), C(m,p)
  integer :: i, j, k
  real :: sum0
  do i = 1, m
    do j = 1, p
      sum0 = 0.0
      do k = 1, n
        sum0 = sum0 + A(i,k) * B(k,j)
      C(i,j) = sum0
      end do
    end do
  end do
end subroutine sub_matmul

!=====================================================================
!======================== MATRICES TRANSPOSE =========================
!=====================================================================

subroutine sub_transpose(a, at)
  implicit none
  real(8), dimension(:,:), intent(in) :: a
  real(8), dimension(size(a, 2), size(a, 1)), intent(out) :: at
  integer :: i, j
  do i = 1, size(a, 1)
    do j = 1, size(a, 2)
      at(j, i) = a(i, j)
    end do
  end do
end subroutine sub_transpose

end program thickplate_q8

!=====================================================================
!======================== END OF THE PROGRAM =========================
!=====================================================================
