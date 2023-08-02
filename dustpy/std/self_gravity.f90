module self_gravity

    implicit none 
    
    contains

    subroutine Omega(r,w)
        ! This subroutine calculates the keplarian angular velocity at a given radius
        !
        ! Parameters:
        ! ------------
        ! r = radial distance from the star (AU)
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
        ! This subroutine calculates the dust opacity at a given temperature (assuming icy grains)
        !
        ! Parameters:
        ! ------------
        ! T      =  midplane gas temperature 
        ! 
        ! Returns:
        ! ------------
        ! kappa  =  opacity
    
        implicit none
    
        integer,  parameter   :: dp = selected_real_kind(15)
        real(dp), intent(in)  :: T
        real(dp), intent(out) :: kappa
    
        real(dp)              :: b = 2
        real(dp)              :: k = 5d-4
    
        kappa = k * T**b
    
    end subroutine 
    
    subroutine f(T,S,out_f)
        ! This subroutine calculates the explicit optical depth given opacity
        !
        ! Parameters:
        ! ------------
        ! T      =  midplane gas temperature 
        ! S      =  midplane gas surface density
        ! 
        ! Returns:
        ! ------------
        ! out_f  =  explicit optical depth 
    
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
        ! This subroutine calculates the temperature where the disk is on the marginally stable
        !
        ! Parameters:
        ! ------------
        ! r  = radial distance from the star (AU)
        ! S  = midplane gas surface density
        ! 
        ! Returns:
        ! ------------
        ! T  = temperature at margin of grav. stability 
        
        use constants
        
        implicit none
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(out) :: T
    
        double precision              :: Q_0 = 1.5
        double precision              :: w
     
        call Omega(r,w)
        T = (2.34 * m_p / k_b) *((pi * G * S * Q_0) / w)**(2.0)
        
    end subroutine T_Q
    
    subroutine FJ_Balance(r,T,S,FJ)
        ! This subroutine calculates the viscous angular momentum flux used to keep the disk in thermal balance
        !
        ! Parameters:
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplane gas temperatur
        ! 
        ! Returns:
        ! ------------
        ! FJ    = viscous angular momentum flux
        
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
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! 
        ! Returns:
        ! ------------
        ! GT    = graviational viscous angular momentum flux
        
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
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplance gas temperature
        ! 
        ! Returns:
        ! ------------
        ! M     = non-graviational angular momentum flux
        
        use constants 
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: M
    
        double precision              :: alpha_m = 0.01
        
        M = 3 * pi * alpha_m * (1.4 * k_B * T / (2.3 * m_p)) * S * r**2.0
    
    end subroutine F_m

    subroutine F_J(r,S,T,FJ)
        ! This subroutine calculates the total angular momentum flux from all sources
        ! Parameters:
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplance gas temperature
        ! 
        ! Returns:
        ! ------------
        ! FJ    = total angular momentum flux

        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: FJ
    
    
        double precision              :: M, GT

        call F_gt(r,S,GT)
        call F_m(r,T,S,M) 

        FJ = max(0.0d0,GT) + M

    end subroutine F_J

    subroutine mass_flux(r,S,T,mass)
        ! This subroutine calculates the accretion rate given a surface density and temperature 
        !
        ! Parameters:
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplance gas temperature
        ! 
        ! Returns:
        ! ------------
        ! mass  =  mass accretion rate   
        
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: mass
    
        double precision              :: FJ, w
    
        call F_J(r,S,T,FJ)
        call Omega(r,w)
    
        !write(*,"(2ES20.5)") GT,M
    
        mass = FJ/(w * r**2.0)
    
    end subroutine mass_flux
    
    subroutine E(r,S,T,En)
        ! This subroutine gives the termal balance relation 
        !
        ! Parameters:
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplance gas temperature
        ! 
        ! Returns:
        ! ------------
        ! En   = mass accretion rate   
        
        implicit none 
    
        double precision, intent(in)  :: r
        double precision, intent(in)  :: S
        double precision, intent(in)  :: T
        double precision, intent(out) :: En
    
        double precision              :: FJ,FJ_B
    
        call F_J(r,S,T,FJ)
        call FJ_Balance(r,T,S,FJ_B)
    
        En = FJ - FJ_B
    
    end subroutine E
    
    subroutine A(r,S,T,diff)
        ! This subroutine gives the difference between the caluclated mass flux and the desired one
        !
        ! Parameters:
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplane gas surface density
        ! T     =  midplance gas temperature
        ! 
        ! Returns:
        ! ------------
        ! diff  = difference in mass accretion rates
        
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
        ! This subroutine calculates good initial guesses for the temp and surface density at a given radius 
        !
        ! Parameters: 
        ! ------------
        ! r     =  radial distance from the star (AU)
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
        ! This subroutine calculates a temperature that keeps the disk in energy balance using the secant method
        !
        ! Parameters: 
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! S     =  midplace gas surface density 
        ! t0    =  first initial temperature guess
        ! t1    =  second initial temperature guess
        !
        ! Returns:
        ! -----------
        ! t2     = final temperature that keeps the disk in energy balance

        use constants
    
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
        ! This subroutine calculates a surface density that keeps the disk in energy balance using the secant method
        !
        ! Parameters: 
        ! ------------
        ! r     =  radial distance from the star (AU)
        ! T     =  midplace gas temperature 
        ! s0    =  first initial surface density guess
        ! s1    =  second initial surface density guess
        !
        ! Returns:
        ! -----------
        ! s2     = final surface density that keeps the disk in energy balance        
    
        implicit none 
    
        double precision, intent(in)     :: r
        double precision, intent(in)     :: T
        double precision, intent(inout)  :: s0
        double precision, intent(inout)  :: s1
        double precision, intent(out)    :: s2
    
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
    
    subroutine init_cond(r,rc,p,S,T,Nr)
        ! This subroutine calculates the initial conditions of the disk by providing the temperature and surface density structure
        ! using the secant method
        !
        ! Parameters: 
        ! ------------
        ! r(Nr)   =  array of radial distances from the star (AU)
        ! rc      =  critical radius 
        ! p       =  exponent 
        ! Nr      =  number of radial grid cells
        !
        ! Returns:
        ! -----------
        ! S(Nr)   = array of midplace gas surface densitys
        ! T(Nr)   = array of midplace gas temperatures

        implicit none 
    
        double precision, intent(in)     :: r(Nr)
        double precision, intent(in)     :: rc
        double precision, intent(in)     :: p
        double precision, intent(out)    :: S(Nr)
        double precision, intent(out)    :: T(Nr)
        integer,          intent(in)     :: Nr
    
        double precision                 :: tau,T_0,T_1,T_2,S_0,S_1,S_2,S_exp
        integer                          :: i
    
        call outer_initial(r(Nr),tau,T_0,S_0)
        T_1 = T_0 + 0.01
        S_1 = S_0 + 0.05
        
        !write(*,"(2ES20.5)") T_0,S_0
    
        do i=Nr,1,-1
            !write(*,"(2ES20.5)") T_0,T_1
            call S_SM(r(i),T_0,S_0,S_1,S_2)

            S_exp = S_2*exp(-(r(i)/rc)**(2+p))

            call T_SM(r(i),S_exp,T_0,T_1,T_2)

            !write(*,"(2ES20.5)") S_exp,T_2
    
            T(i) = T_2
            S(i) = S_exp
    
            T_0 = T_2
            S_0 = S_2
            T_1 = T_0 + 0.01
            S_1 = S_0 + 0.01
    
       end do
    
    end subroutine init_cond
    
    subroutine T_balance(r,S,T_old,T_new,Nr)
        ! This subroutine updates the thermally balanced temperature structure of the disk using the secant method and given a surface density
        !
        ! Parameters: 
        ! ------------
        ! r(Nr)      =  radial distance from the star (AU)
        ! S(Nr)      =  midplace gas surface density
        ! T_old(Nr)  =  previous miplace gas temperature
        ! Nr         =  number of radial grid cells
        !
        ! Returns:
        ! -----------
        ! T_new(Nr)  = thermally balanced temperature structure
        
        implicit none 
    
        double precision, intent(in)     :: r(Nr)
        double precision, intent(in)     :: S(Nr)
        double precision, intent(in)     :: T_old(Nr)
        double precision, intent(out)    :: T_new(Nr)
        integer,          intent(in)     :: Nr            
    
        double precision                 :: T_0,T_1,T_2
        integer                          :: i
        
        T_new(Nr) = 12.0
        T_0 = T_old(Nr)
        !write(*,"(ES20.5)") T_0
        T_1 = T_0 + 0.01
    
        do i=Nr-1,1,-1
            !write(*,"(4ES20.5)") r(i),T_0,T_1,S(i)
            call T_SM(r(i),S(i),T_0,T_1,T_2)
    
            T_0 = T_2
            T_1 = T_0 + 0.01
        end do 
    
    end subroutine T_balance

end module self_gravity