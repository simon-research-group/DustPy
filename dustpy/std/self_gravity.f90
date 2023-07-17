module self_gravity

    implicit none 
    
    contains

    subroutine Omega(r,w)
        ! This subroutine calculates the keplarian angular velocity at a given radius
        !
        ! Parameters:
        ! ------------
        ! r(Nr)  = the array of r values (AU)
        ! Nr     = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! w      = angular velocity 
    
        use constants, only: G,M_sun
        
        implicit none 
    
        integer,  parameter     :: dp = selected_real_kind(15)
        real(dp), intent(in)    :: r
        real(dp), intent(out)   :: w
        ! integer,  intent(in)    :: Nr
    
        w = sqrt(G*M_sun/(r**3))
    
    end subroutine Omega
    
    subroutine opacity(T,kappa)

        implicit none
    
        integer,  parameter   :: dp = selected_real_kind(15)
        real(dp), intent(in)  :: T
        real(dp), intent(out) :: kappa
    
        real(dp)              :: b = 2
        real(dp)              :: k = 2d-4
    
        kappa = k * T**b
    
    end subroutine 
    
    subroutine f(T,S,out_f)
        ! This subroutine calculates the explicit optical depth given opacity
        !
        ! Parameters:
        ! ------------
        ! tau(Nr)    = the array of opacity values 
        ! Nr         = size of the array of tau values
        ! 
        ! Returns:
        ! ------------
        ! out_f(Nr)  = explicit optical depth 
    
        implicit none
    
        integer,  parameter   :: dp = selected_real_kind(15)
        real(dp), intent(in)  :: T
        real(dp), intent(in)  :: S
        real(dp), intent(out) :: out_f
    
        real(dp)              :: kappa, tau
    
        call opacity(T,kappa)
    
        tau = kappa*S
        
        out_f = tau + (1.0/tau)
    
    end subroutine f
    
    subroutine T_Q(r,S,T)
        ! This subroutine calculates the temperature where the disk is on the margin of grav. stability
        !
        ! Parameters:
        ! ------------
        ! r(Nr)  = the array of r values (AU)
        ! Nr     = size of the array of r values
        ! S(Nr)  = the array of surface density values
        ! 
        ! Returns:
        ! ------------
        ! T(Nr)  = temperature at margin of grav. stability 
        
        use constants
        
        implicit none
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(out) :: T
    
        double precision              :: Q_0 = 1.5
        double precision              :: w 
     
        call Omega(r,w)
        T = (2.34 * m_p / k_b) *((pi * G * S * Q_0) / w)**(2.0)
    
        !write(*,"(ES20.5)") T
        
    
    end subroutine T_Q
    
    subroutine FJ_Balance(r,T,S,FJ)
        ! This subroutine calculates the gravitational viscous angular momentum flux 
        !
        ! Parameters:
        ! ------------
        ! r(Nr)   = the array of r values (AU)
        ! S(Nr)   = the array of surface density values
        ! T(Nr)   = the array of temperature values
        ! tau(Nr) = the array of opacity values
        ! Nr      = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! M       = non-graviational angular momentum flux
        
        use constants 
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: T
        double precision, intent(in)  :: S
        double precision, intent(out) :: FJ
    
        double precision              :: T_irr = 12.0
        double precision              :: w
        double precision              :: out_f
        double precision              :: num
        double precision              :: den
        
        call Omega(r,w)
        call f(T,S,out_f)
    
        num = 8 * pi * r**2.0 * sigma_sb * (T**4.0 - T_irr**4.0)
        den = 3 * w * out_f
        
        FJ = num/den
    
    end subroutine FJ_Balance
    
    subroutine F_gt(r,S,GT)
        ! This subroutine calculates the gravitational viscous angular momentum flux 
        !
        ! Parameters:
        ! ------------
        ! r(Nr)   = the array of r values (AU)
        ! S(Nr)   = the array of surface density values
        ! T(Nr)   = the array of temperature values
        ! tau(Nr) = the array of opacity values
        ! Nr      = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! M       = non-graviational angular momentum flux
        
        use constants 
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(out) :: GT
    
        double precision              :: T_irr = 12.0
        double precision              :: T
        double precision              :: w
        double precision              :: out_f
        double precision              :: num
        double precision              :: den
        
        call T_Q(r,S,T)
        call Omega(r,w)
        call f(T,S,out_f)
    
        num = 8.0 * pi * r**2.0 * sigma_sb * (T**4.0 - T_irr**4.0)
        den = 3.0 * w * out_f
        GT = num/den
    
    end subroutine F_gt
    
    
    subroutine F_m(r,T,S,M)
        ! This subroutine calculates the non-gravitational viscous angular momentum flux 
        !
        ! Parameters:
        ! ------------
        ! r(Nr)  = the array of r values (AU)
        ! S(Nr)  = the array of surface density values
        ! T(Nr)  = the array of temperature values
        ! Nr     = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! M      = non-graviational angular momentum flux
        
        use constants 
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: M
    
        double precision              :: alpha_m = 0.01
        
        M = 3 * pi * alpha_m * (k_B * T / (2.34 * m_p)) * S * r**2.0
    
    end subroutine F_m
    
    subroutine mass_flux(r,S,T,mass)
        ! This subroutine calculates the accretion rate given a surface density and temperature 
        !
        ! Parameters:
        ! ------------
        ! r(Nr)   = the array of r values (AU)
        ! S(Nr)   = the array of surface density values
        ! T(Nr)   = the array of temperature values
        ! tau(Nr) = the array of opacity values
        ! Nr      = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! mass    = mass accretion rate   
        
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: mass
    
    
        double precision              :: M
        double precision              :: GT
        double precision              :: w 
    
        call F_gt(r,S,GT)
        call F_m(r,T,S,M)
        call Omega(r,w)
    
        !write(*,"(2ES20.5)") GT,M
    
        mass = (GT + M)/(w * r**2.0)
    
    end subroutine mass_flux
    
    subroutine E(r,S,T,En)
        ! This subroutine calculates the correct temperature that keeps the disk in energy balance
        !
        ! Parameters:
        ! ------------
        ! r(Nr)   = the array of r values (AU)
        ! S(Nr)   = the array of surface density values
        ! T(Nr)   = the array of temperature values
        ! tau(Nr) = the array of opacity values
        ! Nr      = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! mass    = mass accretion rate   
        
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: En
    
    
        double precision              :: M
        double precision              :: GT
        double precision              :: FJ
    
        call F_gt(r,S,GT)
        call F_m(r,T,S,M)
        call FJ_Balance(r,T,S,FJ)
    
        En = (GT + M) - FJ
    
    end subroutine E
    
    subroutine A(r,S,T,diff)
        ! This subroutine gives the difference between the caluclated mass flux and the desired one
        !
        ! Parameters:
        ! ------------
        ! r(Nr)   = the array of r values (AU)
        ! S(Nr)   = the array of surface density values
        ! T(Nr)   = the array of temperature values
        ! tau(Nr) = the array of opacity values
        ! Nr      = size of the array of r values
        ! 
        ! Returns:
        ! ------------
        ! mass    = mass accretion rate   
        
        use constants, only: M_sun,year
    
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: diff
    
        double precision              :: mass
    
        call mass_flux(r,S,T,mass)
    
        diff = mass - (7d-5 * M_sun /year)
    
    end subroutine A
    
    subroutine outer_initial(r,tau,T,S)
        ! This subroutine calculates good initial guesses for the temp and surface density at 300 AU
        !
        ! Parameters: 
        ! ------------
        ! r(Nr)   = array of r values (AU)
        ! Nr      = size of the r array
        !
        ! Returns:
        ! -----------
        ! tau     = initial guess for opacity 
        ! T_0     = initial guess for temp
        ! S_0     = initial guess for the surface density 
    
        use constants
    
        implicit none
    
        double precision, intent(in)  :: r
        double precision, intent(out) :: tau
        double precision, intent(out) :: T
        double precision, intent(out) :: S
        
        double precision              :: tau_vals(12)
        double precision              :: T_0(12)
        double precision              :: S_0(12)
        double precision              :: w
        double precision              :: m_dots(12)
        double precision              :: T_irr = 12.0
        double precision              :: Q_0   = 1.5
        integer                       :: index,i  
    
    
        do i=1,12
            tau_vals(i) = i*0.1
        end do
        
        call Omega(r,w)
    
        T_0(:) = (T_irr**4 + (3*(tau_vals(:)+1.0/tau_vals(:))*(7d-5*M_sun/year)*w**2)/(8*pi*sigma_sb))**(1.0/4.0)
    
        S_0(:) = w/(pi * G * Q_0) * (k_B * T_0(:)/(2.34 * m_p))**(1.0/2.0)
        
        do i=1,12
            call mass_flux(r,S_0(i),T_0(i),m_dots(i))
        end do
    
        index = minloc(abs((7d-5 * M_sun /year) - m_dots), dim = 1)
    
        tau = tau_vals(index)
        T = T_0(index)
        S = S_0(index)
    
    end subroutine outer_initial
    
    subroutine T_SM(r,S,t0,t1,t2)
    
        implicit none 
    
        double precision, intent(in)     :: r
        double precision, intent(in)     :: S
        double precision, intent(inout)  :: t0
        double precision, intent(inout)  :: t1
        double precision, intent(out)    :: t2
    
        double precision  :: tol = 1d-12
        double precision  :: E0,E1
    
        call E(r,S,t1,E1)
        do while (ABS(E1)>=tol .and. ABS(t1-t0)>=tol)
            call E(r,S,t0,E0)
            call E(r,S,t1,E1)
    
            t2 = t1 - (E1)*(t1 - t0)/(E1 - E0)
            t0 = t1
            t1 = t2
        
        end do 
    
    end subroutine T_SM
    
    subroutine S_SM(r,T,s0,s1,s2)
    
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: T
        double precision, intent(inout)  :: s0
        double precision, intent(inout)  :: s1
        double precision, intent(out) :: s2
    
        double precision  :: tol = 1d-12
        double precision  :: A0,A1
    
        call A(r,s1,T,A1)
        do while (ABS(A1)>=tol .and. ABS(s1-s0)>=tol)
            call A(r,s0,T,A0)
            call A(r,s1,T,A1)
    
            s2 = s1 - (A1)*(s1 - s0)/(A1 - A0)
    
            s0 = s1
            s1 = s2
        end do 
    
    end subroutine S_SM
    
    subroutine init_cond(r,S,T,Nr)
        implicit none 
    
        double precision, intent(in)     :: r(Nr)
        double precision, intent(out)    :: S(Nr)
        double precision, intent(out)    :: T(Nr)
        integer,          intent(in)     :: Nr
    
        double precision                 :: tau,T_0,T_1,T_2,S_0,S_1,S_2
        integer                          :: i
    
        call outer_initial(r(Nr),tau,T_0,S_0)
        T_1 = T_0 + 0.01
        S_1 = S_0 + 0.05
    
        do i=Nr,1,-1
            call S_SM(r(i),T_0,S_0,S_1,S_2)
            call T_SM(r(i),S_2,T_0,T_1,T_2)
    
            T(i) = T_2
            S(i) = S_2
    
            T_0 = T_2
            S_0 = S_2
            T_1 = T_0 + 0.01
            S_1 = S_0 + 0.05
    
       end do
    
    end subroutine init_cond
    
    subroutine T_balance(r,S,T_old,T_new,Nr)
    
        implicit none 
    
        double precision, intent(in)     :: r(Nr)
        double precision, intent(in)     :: S(Nr)
        double precision, intent(in)     :: T_old(Nr)
        double precision, intent(out)    :: T_new(Nr)
        integer,          intent(in)     :: Nr            
    
        double precision                 :: T_0,T_1,T_2
        integer                          :: i
        
        T_0 = T_old(Nr)
        T_1 = T_0 + 0.01
    
        do i=Nr,1,-1
            call T_SM(r(i),S(i),T_0,T_1,T_2)
            T_new(i) = T_2
    
            T_0 = T_2
            T_1 = T_0 + 0.01
        end do 
    
    end subroutine T_balance

end module self_gravity