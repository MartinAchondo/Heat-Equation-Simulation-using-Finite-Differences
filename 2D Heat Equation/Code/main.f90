
! Programa para resolver la ecuacion de calor
program main
    implicit none
    integer :: Nx,Ny,tout,nt
    real(8) :: time_f,d_time,dx,dy,Lx,Ly,q
    real(8), allocatable, dimension(:,:) :: Temperature_f,Temperature_i,qn,qx,qy
    real(8), allocatable, dimension(:) :: Res
    real(8) :: alpha,k,sigma,theta,cp,rho,sigma2
    character(25) :: file_T_eu,file_T_cn,file_qxeu,file_qyeu,file_q_cn,file_q_eu,file_qxcn,file_qycn,file_reseu,file_rescn
    logical :: euler,crank,flag_results
    
    ! Codigo exporta los siguientes archivos con resultados

    file_T_eu = 'Results_T_eu.txt'
    file_T_cn = 'Results_T_cn.txt'
    file_q_eu = 'Results_q_eu.txt'
    file_q_cn = 'Results_q_cn.txt'
    file_qxeu = 'Results_qx_eu.txt'
    file_qyeu = 'Results_qy_eu.txt'
    file_qxcn = 'Results_qx_cn.txt'
    file_qycn = 'Results_qy_cn.txt'
    file_reseu = 'Residual_eu.txt'
    file_rescn = 'Residual_cn.txt'

    k = 100.0d0
    cp = 421.0d0
    rho = 8862.0d0
    alpha = k/(cp*rho)
    q = 1.0d4

    time_f = 60000
    d_time = 1
    tout = 1000
    nt = floor(time_f/d_time)+1
    dx = 0.1d0
    dy = dx
    Lx = 2d0
    Ly = 1d0

    sigma = alpha*d_time/(dx**2)
    sigma2 = 2*sigma
    theta = -2*dx/k

    Nx = floor(Lx/dx) + 1
    Ny = floor(Ly/dy) + 1
    allocate (Temperature_i(Ny,Nx))
    allocate (Temperature_f(Ny,Nx))
    allocate (qn(Ny,Nx))
    allocate (qx(Ny,Nx))
    allocate (qy(Ny,Nx))
    allocate (Res(nt))

    Temperature_i = 0
    
    write(*,*)
    write(*,*) "Codigo resolucion ecuacion de calor placa"
    write(*,*)
    write(*,*) "Propiedades (k,alpha,rho,cp)"
    write(*,"(F7.1,ES12.4,F8.1,F7.1)") k,alpha,rho,cp
    write(*,*) "Numero de Fourier: ",sigma2
    write(*,*) "Nodos (Nx,Ny):",Nx,Ny
    write(*,*) "Espaciamientos (dt,dx=dy)"
    write(*,"(2F10.5)") d_time,dx

    write(*,*)
    write(*,*)

    call system('mkdir eu_rtxt')
    call system('mkdir cn_rtxt')
    call system('mkdir results')

    write(*,*)
   
    ! Condiciones para correr esquemas
    euler = .true.
    crank = .true.

    ! Condicion para imprimir los resultados cada tout
    flag_results = .false.

    if (euler) then
        write(*,*)
        write(*,*)
        write(*,*) "------- EULER EXPLICITO --------"
        write(*,*)

        call Euler_Explicito(Nx,Ny,Temperature_i,Temperature_f,d_time,time_f,sigma,q,theta,tout,flag_results,Res,nt)
        call WriteFile(Nx,Ny,Temperature_f,file_T_eu,.true.,'results/')
        call Calculate_Heat(Nx,Ny,Temperature_f,qn,theta,q,qx,qy)
        call WriteFile(Nx,Ny,qn,file_q_eu,.true.,'results/')
        call WriteFile(Nx,Ny,qx,file_qxeu,.true.,'results/')
        call WriteFile(Nx,Ny,qy,file_qyeu,.true.,'results/')
        call WriteFile2(nt,Res,file_reseu,.true.,'results/')

    end if

    if (crank) then
        write(*,*)
        write(*,*)
        write(*,*) "------- CRANK NICOLSON --------"
        write(*,*)
        
        call Crank_Nicolson(Nx,Ny,Temperature_i,Temperature_f,d_time,time_f,sigma,q,theta,tout,flag_results,Res,nt)
        call WriteFile(Nx,Ny,Temperature_f,file_T_cn,.true.,'results/')
        call Calculate_Heat(Nx,Ny,Temperature_f,qn,theta,q,qx,qy)
        call WriteFile(Nx,Ny,qn,file_q_cn,.true.,'results/')
        call WriteFile(Nx,Ny,qn,file_qxcn,.true.,'results/')
        call WriteFile(Nx,Ny,qn,file_qycn,.true.,'results/')
        call WriteFile2(nt,Res,file_rescn,.true.,'results/')
    end if

    write(*,*)
    write(*,*)
    write(*,*)

end program main


!============ Euler Explicito ===================!

subroutine Euler_Explicito(Nx,Ny,Ti,Tf,dt,timef,sigma,q,theta,tout,flag,Res,ntt)
    implicit none
    integer, intent(in) :: Nx,Ny,tout,ntt
    real(8), intent(in) :: dt,timef,sigma,Ti(Ny,Nx),q,theta
    real(8), intent(out) :: Tf(Ny,Nx),Res(ntt)
    logical, intent(in) :: flag
    integer :: i,j,nt,kk,ll,stat
    real(8) :: T_old(Ny,Nx),T(Ny,Nx),time,T_est(Ny,Nx),q2,Calc_Residual
    character(7) :: x1

    nt = 0

    T_old(:,:) = Ti(:,:)
    time = 0
    do while (time<timef) 
        do i=2,Nx-1
            do j=2,Ny-1
                T(j,i)=T_old(j,i)+sigma*(T_old(j,i-1)-2*T_old(j,i)+T_old(j,i+1))+sigma*(T_old(j-1,i)-2*T_old(j,i)+T_old(j+1,i))
            end do
        end do
        do i=2,Nx-1
            if (i.le.Ny) then
                q2 = q
            else
                q2 = 0d0
            end if 
            j=1
            T(j,i)=T_old(j,i)+sigma*(T_old(j,i-1)-2*T_old(j,i)+T_old(j,i+1))+sigma*(-2*T_old(j,i)+2*T_old(j+1,i)-q2*theta)
            j=Ny
            T(j,i)=T_old(j,i)+sigma*(T_old(j,i-1)-2*T_old(j,i)+T_old(j,i+1))+sigma*(-2*T_old(j,i)+2*T_old(j-1,i)-q2*theta)
        end do

        T(:,1) = 0
        T(:,Nx) = 0

        time = time + dt

        ! Criterio estacionario
        T_est = T - T_old
        if (maxval(abs(T_est)).lt.0.0001) then
            write(*,*) 'Estado estacionario, tiempo:', time
            write(*,*) 'Max Temp:',maxval(T)
            exit
        end if

        ! Criterio si imprimir cada tout
        if (flag) then
            if (mod(nt,tout)==0 .or. nt==0) then 
                write (x1,'(I7.7)') floor(time)
                open(unit=20,file='eu_rtxt/t_eu'//trim(x1)//'.txt',iostat=stat,action='write')
                if(stat/=0) then
                    write(*,*)'Error al abrir el archivo con iostat',stat
                end  if
                do kk=1,ny
                    write(20,*)(T(kk,ll), ll=1,nx)
                end do
                close(unit=20,iostat=stat)
            end if
        end if

        ! Residual
        Res(nt+1) = Calc_Residual(Nx,Ny,T_old,T,sigma,dt)

        nt = nt + 1
        T_old(:,:) = T(:,:)
    end do
    Tf(:,:) = T_old(:,:)


end subroutine Euler_Explicito


!===== Crank Nicolson ==============!

subroutine Crank_Nicolson(Nx,Ny,Ti,Tf,dt,timef,sigma,q,theta,tout,flag,Res,ntt)
    implicit none
    integer, intent(in) :: Nx,Ny,tout,ntt
    real(8), intent(in) :: dt,timef,sigma,Ti(Ny,Nx),q,theta
    real(8), intent(out) :: Tf(Ny,Nx),Res(ntt)
    logical, intent(in) :: flag
    integer :: Nv,nt,kk,ll,stat
    character(7) :: x1
    real(8) :: T_old(Ny,Nx),T(Ny,Nx),time,Calc_Residual,T_est(Ny,Nx)

    real(8), allocatable, dimension(:) :: b_vec,T_vec
    real(8), allocatable, dimension(:,:) :: A,A2

    Nv = (Nx-2)*(Ny)
    allocate(b_vec(Nv))
    allocate(T_vec(Nv))
    allocate(A(Nv,Nv))
    allocate(A2(Nv,Nv))

    ! Crear matriz
    call create_mat(A,sigma,Nv,Nx,Ny)
    write(*,*) 'Matrix created CN'

    T_old(:,:) = Ti(:,:)
    time = 0
    nt = 0
    do while (time<timef) 
        ! Create b_system
        call b_system(T_old,b_vec,Nx,Ny,Nv,sigma,q,theta)
        ! resolver sistema
        call Gseid(A,b_vec,Nv,T_vec,300,1.0d-8,0.5d0)
        ! vec to mat T
        call vec_to_mat(T,T_vec,Nx,Ny,Nv)

        time = time + dt
        
        ! Criterio de estacionario
        T_est = T - T_old
        if (maxval(abs(T_est)).lt.0.0001) then
            write(*,*) 'Estado estacionario, tiempo:', time
            write(*,*) 'Max Temp:',maxval(T)
            exit
        end if

        ! Criterio de imprimir cada tout
        if (flag) then
            if (mod(nt,tout)==0 .or. nt==0) then 
                write (x1,'(I7.7)') floor(time)
                open(unit=20,file='cn_rtxt/t_eu'//trim(x1)//'.txt',iostat=stat,action='write')
                if(stat/=0) then
                    write(*,*)'Error al abrir el archivo con iostat',stat
                end  if
                do kk=1,ny
                    write(20,*)(T(kk,ll), ll=1,nx)
                end do
                close(unit=20,iostat=stat)
            end if
        end if

        ! Residual
        Res(nt+1) = Calc_Residual(Nx,Ny,T_old,T,sigma,dt)

        nt = nt + 1
        T_old(:,:) = T(:,:)
    end do
    Tf(:,:) = T_old(:,:)

end subroutine Crank_Nicolson


! Funcion que calcula el residual en un paso
real(8) function Calc_Residual(Nx,Ny,Told,Tnew,sigma,dt)
    implicit none
    integer, intent(in) :: Nx,Ny
    real(8), intent(in) :: Told(Ny,Nx),Tnew(Ny,Nx),sigma,dt
    integer :: i,j
    real(8) :: Res(Ny-2,Nx-2),rxy

    do j=2,Ny-1
        do i=2,Nx-1
            Res(j-1,i-1)=(Tnew(j,i)-Told(j,i))/dt-(sigma/dt)*(Tnew(j,i-1)-4*Tnew(j,i)+Tnew(j,i+1)+Tnew(j-1,i)+Tnew(j+1,i))
        end do
    end do

    rxy = 0.0d0
    do j=1,Ny-2
        do i=1,Nx-2
            rxy = rxy + (1.0d0/((Nx-2)*(Ny-2)))*Res(j,i)**2
        end do
    end do

    Calc_Residual = sqrt(rxy)
    
end function Calc_Residual


! Creacion de matriz para CN
subroutine create_mat(A,sigma,Nv,Nx,Ny)
    implicit none
    integer, intent(in) :: Nv,Nx,Ny
    real(8), intent(in) :: sigma
    real(8), intent(out) :: A(Nv,Nv)
    integer :: i,j,k,nxx,nyy
    real(8) :: mu

    nxx = Nx-2
    nyy = Ny
    A = 0
    mu = sigma/2.0d0

    do j=1,nyy
        do i=1,nxx
            k = i + (j-1)*nxx
            if (j==1) then
                if (i==1) then
                    A(k,k) = 1 + 4*mu
                    A(k,k+1) = -mu
                elseif (i==nxx) then
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                else
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                    A(k,k+1) = -mu
                end if
                A(k,k+nxx) = -2*mu
            elseif (j==nyy) then
                if (i==1) then
                    A(k,k) = 1 + 4*mu
                    A(k,k+1) = -mu
                elseif (i==nxx) then
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                else
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                    A(k,k+1) = -mu
                end if
                A(k,k-nxx) = -2*mu
            else
                if (i==1) then
                    A(k,k) = 1 + 4*mu
                    A(k,k+1) = -mu
                elseif (i==nxx) then
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                else
                    A(k,k) = 1 + 4*mu
                    A(k,k-1) = -mu
                    A(k,k+1) = -mu
                end if
                A(k,k-nxx) = -mu
                A(k,k+nxx) = -mu
            end if
        end do
    end do

end subroutine create_mat


! Creacion de vector b para Ax=b en CN
subroutine b_system(T_mat,b_vec,Nx,Ny,Nv,sigma,q,theta)
    implicit none
    integer, intent(in) :: Nx,Ny,Nv
    real(8), intent(in) :: T_mat(Ny,Nx),sigma,q,theta
    real(8), intent(out) :: b_vec(Nv)
    integer :: i,j,k
    real(8) :: mu,q2

    mu = sigma/2.0d0
    k = 0
    do j=1,Ny
        do i=2,Nx-1
            k = k + 1
            if (i.le.Ny) then
                q2 = q
            else
                q2 = 0
            end if
            if (j==1) then
                b_vec(k)=mu*T_mat(j,i+1)+(1.0d0-4.0d0*mu)*T_mat(j,i)+mu*T_mat(j,i-1)+2.0d0*mu*T_mat(j+1,i)-2*mu*q2*theta
            else if (j==Ny) then
                b_vec(k)=mu*T_mat(j,i+1)+(1.0d0-4.0d0*mu)*T_mat(j,i)+mu*T_mat(j,i-1) +2.0d0*mu*T_mat(j-1,i)-2*mu*q2*theta
            else
                b_vec(k) = mu*T_mat(j,i+1) + (1.0d0-4.0d0*mu)*T_mat(j,i) + mu*T_mat(j,i-1) + mu*T_mat(j+1,i)+mu*T_mat(j-1,i) 
            end if
        end do
    end do

end subroutine b_system

! Paso del vector solucion a matriz de temperaturas en CN
subroutine vec_to_mat(T_mat,T_vec,Nx,Ny,Nv)
    implicit none
    integer, intent(in) :: Nx,Ny,Nv
    real(8), intent(in) :: T_vec(Nv)
    real(8), intent(out) :: T_mat(Ny,Nx)
    integer :: i,j,k

    T_mat = 0.0d0
    k = 0
    do j=1,Ny
        do i=2,Nx-1
            k = k + 1
            T_mat(j,i) = T_vec(k)
        end do
    end do

end subroutine vec_to_mat


!============ Calcular Flujo de Calor ============!
! Uitiliza diferencias centradas en nodos interiores y adelantadas o atrasadas en bordes.

subroutine Calculate_Heat(Nx,Ny,T,qn,theta,qw,qxx,qyy)
    implicit none
    integer, intent(in) :: Nx,Ny
    real(8), intent(in) :: theta,T(Ny,Nx),qw
    real(8), intent(out) :: qn(Ny,Nx),qxx(Ny,Nx),qyy(Ny,Nx)
    integer :: i,j
    real(8) :: qx,qy,qq(Ny,Nx)

    qq = 0
    do i=2,Nx-1
        do j=2,Ny-1
            qx = (1d0/theta)*(T(j,i+1)-T(j,i-1))
            qy = (1d0/theta)*(T(j+1,i)-T(j-1,i))
            qq(j,i) = sqrt(qx**2+qy**2)
            qxx(j,i) = qx
            qyy(j,i) = qy
        end do
    end do
    do i=2,Nx-1
        if (i.lt.Ny) then
            qx = (1d0/theta)*(T(1,i+1)-T(1,i-1))
            qy = qw
        else
            qx = (1d0/theta)*(T(1,i+1)-T(1,i-1))
            qy = 0
        end if
        qq(1,i) = sqrt(qx**2+qy**2)
        qxx(1,i) = qx
        qyy(1,i) = qy
    end do
    do i=2,Nx-1
        if (i.lt.Ny) then
            qx = (1d0/theta)*(T(Ny,i+1)-T(Ny,i-1))
            qy = -qw
        else
            qx = (1d0/theta)*(T(Ny,i+1)-T(Ny,i-1))
            qy = 0
        end if
        qq(Ny,i) = sqrt(qx**2+qy**2)
        qxx(Ny,i) = qx
        qyy(Ny,i) = qy
    end do
    do j=2,Ny-1
        qx = (1d0/theta)*(-3*T(j,1)+4*T(j,2)-T(j,3))
        qy = (1d0/theta)*(T(j+1,1)-T(j-1,1))
        qq(j,1) = sqrt(qx**2+qy**2)
        qxx(j,1) = qx
        qyy(j,1) = qy
    end do
    do j=2,Ny-1
        qx = (1d0/theta)*(3*T(j,Nx)-4*T(j,Nx-1)+T(j,Nx-2))
        qy = (1d0/theta)*(T(j+1,Nx)-T(j-1,Nx))
        qq(j,Nx) = sqrt(qx**2+qy**2)
        qxx(j,1) = qx
        qyy(j,1) = qy
    end do
    qn(:,:) = qq(:,:)
end subroutine Calculate_Heat

! Metodo de Gauss Seidel con relajacion
subroutine Gseid(a_mat,b2,n,x,imax,es,lambda)
    implicit none
    integer, intent(in) :: n,imax
    real(8), intent(in) :: a_mat(n,n),b2(n),es,lambda
    real(8), intent(out) :: x(n)
    integer :: iter,centinela,i,j
    real(8) :: a(n,n),b(n),dummy,sum,old,ea

    a = a_mat
    b = b2

    do i=1,n
        dummy = a(i,i)
        do j=1,n
            a(i,j) = a(i,j)/dummy
        end do
        b(i) = b(i)/dummy
    end do
    do i=1,n
        sum = b(i)
        do j=1,n
            if(i.ne.j) then
                sum = sum - a(i,j)*x(j)
            end if
        end do
        x(i) = sum
    end do
    iter = 1
    do while(.true.)
        centinela = 1
        do i=1,n
            old = x(i)
            sum = b(i)
            do j=1,n
                if(i.ne.j) then
                    sum = sum - a(i,j)*x(j)
                end if
            end do
            x(i) = lambda*sum + (1-lambda)*old
            if (centinela == 1 .and. x(i) .ne.0 ) then
                ea = abs((x(i)-old)/x(i))*100
                if (ea.gt.es) then
                    centinela = 0
                end if
            end if
        end do
        iter = iter + 1
        if (centinela==1 .or. iter.ge.imax) then
            exit
        end if
    end do

end subroutine Gseid


! Para escribir archivo salida
subroutine WriteFile(nx,ny,T,file,flag,dir)
    implicit none
    integer, intent(in) :: nx,ny
    character(25), intent(in) :: file
    character(8), intent(in) :: dir
    logical, intent(in) :: flag
    integer :: i,j,stat
    real(kind=8), intent(in) :: T(ny,nx)

    open(unit=20,file=dir//file,iostat=stat,action='write')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do j=1,ny
            write(20,*)(T(j,i), i=1,nx)
        end do
    close(unit=20,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    if (flag) then
        write(*,*) "Archivo Guardado: "//file
    end if

end subroutine WriteFile

! Para escribir archivo de salida
subroutine WriteFile2(nt,V,file,flag,dir)
    implicit none
    integer, intent(in) :: nt
    character(25), intent(in) :: file
    character(8), intent(in) :: dir
    logical, intent(in) :: flag
    integer :: j,stat
    real(kind=8), intent(in) :: V(nt)

    open(unit=20,file=dir//file,iostat=stat,action='write')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do j=1,nt
            write(20,*)(V(j))
        end do
    close(unit=20,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    if (flag) then
        write(*,*) "Archivo Guardado: "//file
    end if

end subroutine WriteFile2


