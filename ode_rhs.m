function [ output_args ] = ode_rhs( x, z )

global P0 A mdot_ox mdot_l0 rho_l FO_st hfg ii D0

D = sqrt(z(1))
Tg = z(2)
vd = z(3)
mdot_g = z(4)
phi_g = z(5)
x

    if ii == 1
       D0 = D;
    end
    
    ii = ii + 1;

    if D > 0 && isreal(D) == 1
        gas = GRI30('Multi');
        nsp = nSpecies(gas);
        iLOX = speciesIndex(gas,'O2');
        iCH4 = speciesIndex(gas,'CH4');

        FO = phi_g * FO_st;
        y = zeros(nsp,1);
        y(iCH4,1) = FO / (1 + FO);
        y(iLOX,1) = 1 / (1 + FO);

        set(gas,'Temperature',Tg,'Pressure',P0,'Y',y);
        rho_g = density(gas);
        %equilibrate(gas,'HP');
        mu_g = viscosity(gas);

        vg = mdot_g / A / rho_g;
        v_rel = vg - vd;
        Re = rho_g * v_rel * D / mu_g;
        Cd = 24/Re + 6/(1+sqrt(Re)) + 0.4;

        K = K_methane(Tg,P0,phi_g);

        %dml(ii) = -1.5 * mdot_l0 * D * K / (D0^3 * vd);
        %dphi(ii) = -dml(ii) / (FO_st * mdot_ox);
        dml = -1.5 * mdot_l0 * D * K / (D0^3 * vd);
        dphi = -dml/(FO_st * mdot_ox);
        
        h1 = dhg_dphi(Tg,P0,phi_g);
        h2 = dhg_dT(Tg,P0,phi_g);

        output_args = [-K/vd;
            (hfg/mdot_g*dml-h1*dphi)/h2;
            3*Cd*rho_g*v_rel*abs(v_rel) / (4*rho_l*D*vd);
            -dml;
            dphi];
    else
        output_args = [0; 0; 0; 0; 0];
    end
end

