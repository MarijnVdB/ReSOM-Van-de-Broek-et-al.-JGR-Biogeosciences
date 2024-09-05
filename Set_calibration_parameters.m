%% This script is used to set the calibration parameters, depending on the thermal adaptation scenario

if thermalAdapt == 0
    
    % The parameters that need to be calibrated are defined with empty values
    res_turnover_rate_Rstrat = []; % 1
    micr_death_rate_Rstrat = []; % 2
    max_poly_degradation_Rstrat = []; % 3
    max_enz_prod_Rstrat = []; % 4
    max_doc_uptake_Rstrat = []; % 5
    max_monomer_ads_organic = []; % 6

    micr_maint_rate_Kstrat = []; % 7
    micr_death_rate_Kstrat = []; % 8
    max_micr_growth_Kstrat = []; % 9

    res_turnover_rate = []; % 10
    micr_death_rate = []; % 11
    max_poly_degradation = []; % 12
    max_doc_uptake = []; % 13
    max_monomer_ads_soil = []; % 14

    % Parameters for the temperature response function
    par_MMRT_1 = []; % 15
    Topt = []; % 16

    % An array with these empty values is constructed
    parameters = [res_turnover_rate_Rstrat; micr_death_rate_Rstrat; max_poly_degradation_Rstrat; ...
        max_enz_prod_Rstrat; max_doc_uptake_Rstrat; max_monomer_ads_organic; micr_maint_rate_Kstrat; micr_death_rate_Kstrat; max_micr_growth_Kstrat; ...
        res_turnover_rate; micr_death_rate; max_poly_degradation; max_doc_uptake; max_monomer_ads_soil; ...
        par_MMRT_1; Topt];

    % The upper and lower bounds of the parameter values are defined,
    % in the same order in which the parameters are defined above
    lb = [0.1; 0.01;  1; 0.0005;  50; 1e-4; 0.005; 0.01; 0.01;  1; 0.005;   1; 200; 1e-4; -6; 273.15+25];
    ub = [  1;  0.1; 50;   0.01; 500; 1e-1;   0.1;  0.1;    2; 50;   0.1; 150; 600; 1e-1; -1; 273.15+50];

    % An initial population of parameter values
%     initialPop = [0.3033, 0.0616, 37.1554, 0.0018, 273.3270, 0.004, 0.0547, 0.0366, 0.0277, 5.6649, 0.0413, 93.0579, 475.0144, 0.01, -2.17, 306.2];
    initialPop = [0.973533503584554,0.0897860739768491,13.7151549875821,0.00164816592543630,382.703487663186,0.00191593561030590,0.0129679832868755,0.0725322960197611,0.713776363964896,40.4480981993862,0.0726794065500046,66.4027755044967,412.568986531371,0.00313586346374902,-4.23719949375680,301.390298380053];

    % The number of parameters is calculated and stored
    numberOfParametersToCalibrate = numel(lb);
    
    % If theta_op is calibrated, this parameter is added
    if calib_moisture == 1
        
        theta_op = []; % 17
        
        parameters = [parameters; theta_op];
        lb = [lb; 0.1];
        ub = [ub; 0.9];
        initialPop = [0.922003625351503,0.0661654079356321,12.2151549875821,0.00218722337480572,449.714055632997,0.00191593561030590,0.0381345568084556,0.0412205671258852,0.963776363964896,40.4480981993862,0.0203493945195814,32.8511324026577,593.383007151986,0.00840727408567065,-2.17006373351896,311.811991200732,0.315237451802814];
%         clear('initialPop')
        numberOfParametersToCalibrate = numel(lb);
        
    end
    
elseif thermalAdapt == 1
    
    if optimumDriven == 1

        res_turnover_rate_Rstrat = []; % 1
        micr_death_rate_Rstrat = []; % 2
        max_poly_degradation_Rstrat = [];  % 3
        max_enz_prod_Rstrat = []; % 4
        max_doc_uptake_Rstrat = []; % 5
        max_monomer_ads_organic = []; % 6

        micr_maint_rate_Kstrat = []; % 7
        micr_death_rate_Kstrat = []; % 8
        max_micr_growth_Kstrat = []; % 9

        res_turnover_rate = []; % 10
        micr_death_rate = []; % 11
        max_poly_degradation = []; % 12
        max_doc_uptake = []; % 13
        max_monomer_ads_soil = []; % 14

        % Parameters for the temperature response function
        par_MMRT_1 = []; % 15
        alphaT = []; % 16
        betaT = []; % 17
        nYears_thermalAdapt = []; % 18

        % An array with these empty values is constructed
        parameters = [res_turnover_rate_Rstrat; micr_death_rate_Rstrat; max_poly_degradation_Rstrat; ...
            max_enz_prod_Rstrat; max_doc_uptake_Rstrat; max_monomer_ads_organic; micr_maint_rate_Kstrat; micr_death_rate_Kstrat; max_micr_growth_Kstrat; ...
            res_turnover_rate; micr_death_rate; max_poly_degradation; max_doc_uptake; max_monomer_ads_soil; ...
            par_MMRT_1; alphaT; betaT; nYears_thermalAdapt];

        % The upper and lower bounds of the parameter values are defined,
        % in the same order in which the parameters are defined above
        lb = [0.1; 0.01;  1; 0.0005;  50; 1e-4; 0.005; 0.01; 0.01;  1; 0.005;   1; 200; 1e-4; -6; 0; 10; 1/6];
        ub = [  1;  0.1; 50;   0.01; 500; 1e-1;   0.1;  0.1;    2; 50;   0.1; 100; 600; 1e-1; -1; 1; 50;   3];

        % An initial population of parameter values
        initialPop = [0.672003625351503,0.0802204691636239,11.7151549875821,0.00189189165084145,402.119673309675,0.00191593561030590,0.0604356515548861,0.0376971023827847,0.527192017089541,39.4480981993862,0.0328699537307972,36.3461181956435,589.383007151986,0.00840727408567065,-5.28965646622816,0.646346065444023,20.2008265155231,0.946390063061803];

        % The number of parameters is calculated and stored
        numberOfParametersToCalibrate = numel(lb);
        
        % If theta_op is calibrated, this parameter is added
        if calib_moisture == 1

            theta_op = []; % 19

            parameters = [parameters; theta_op];
            lb = [lb; 0.1];
            ub = [ub; 0.9];
            initialPop = [0.988428347772265,0.0897860739768491,38.3216392292561,0.000661337914203982,459.310634900607,0.00191593561030590,0.00803264980287795,0.0591898413774663,0.287096461161694,6.55673828245017,0.00500000000000000,50.7993902687911,590.383007151986,0.0157651079656851,-4.58400683079527,0.265301726053799,26.2238080128241,1.01895909020539,0.470635903263446];
%             clear('initialPop')
            numberOfParametersToCalibrate = numel(lb);

        end
        
    elseif enzymeRigidity == 1
        
        res_turnover_rate_Rstrat = []; % 1
        micr_death_rate_Rstrat = []; % 2
        max_poly_degradation_Rstrat = []; % 3
        max_enz_prod_Rstrat = []; % 4
        max_doc_uptake_Rstrat = []; % 5
        max_monomer_ads_organic = []; % 6

        micr_maint_rate_Kstrat = []; % 7
        micr_death_rate_Kstrat = []; % 8
        max_micr_growth_Kstrat = []; % 9

        res_turnover_rate = []; % 10
        micr_death_rate = []; % 11
        max_poly_degradation = []; % 12
        max_doc_uptake = []; % 13
        max_monomer_ads_soil = []; %14

        alphaT = []; % 15
        betaT = []; % 16
        alphaC = []; % 17
        betaC = []; % 18
        nYears_thermalAdapt = []; % 19

        % An array with these empty values is constructed
        parameters = [res_turnover_rate_Rstrat; micr_death_rate_Rstrat; max_poly_degradation_Rstrat; ...
            max_enz_prod_Rstrat; max_doc_uptake_Rstrat; max_monomer_ads_organic; micr_maint_rate_Kstrat; micr_death_rate_Kstrat; max_micr_growth_Kstrat; ...
            res_turnover_rate; micr_death_rate; max_poly_degradation; max_doc_uptake; max_monomer_ads_soil; ...
            alphaT; betaT; alphaC; betaC; nYears_thermalAdapt];

        % The upper and lower bounds of the parameter values are defined,
        % in the same order in which the parameters are defined above
        lb = [0.1; 0.01;  1; 0.0005;  50; 1e-4; 0.005; 0.01; 0.01;  1; 0.005;   1; 200; 1e-4; 0; 10; 0; -10; 1/6];
        ub = [  1;  0.1; 50;   0.01; 500; 1e-1;   0.1;  0.1;    2; 50;   0.1; 100; 600; 1e-1; 1; 50; 1;  -1;   3];

        % Some parameters have to be larger than others
        A = [zeros(1,16), 12, 1, 0];
        b = 0;

        % An initial population of parameter values
%         initialPop = [0.3456, 0.0466, 4.9036, 0.0060, 91.2815, 0.002, 0.0140, 0.0394, 0.8220, 2.3263, 0.0227, 14.0622, 208.6907, 0.002, 0.8607, 288.4560, 0.3587, -12.7752, 1];
        initialPop = [0.523484955408273,0.0503672168124505,2.67058555131025,0.00599759597654422,92.6143386716884,0.00140253021123305,0.0141179180698909,0.0428105196300750,0.776856159272648,2.65660648872450,0.0271590877124212,10.6928253943859,210.223436970382,0.00348246783540101,0.856852169631644,15.5206901368021,0.334738047883527,-9.61658499917097,0.999140950870731];
        
        % The number of parameters is calculated and stored
        numberOfParametersToCalibrate = numel(lb);
        
        % If theta_op is calibrated, this parameter is added
        if calib_moisture == 1

            theta_op = []; % 20

            parameters = [parameters; theta_op];
            lb = [lb; 0.1];
            ub = [ub; 0.9];
            initialPop = [0.694229651549651,0.0985757242000947,2.62272656456839,0.00980626620621209,354.434059819392,0.00177077850494300,0.00800291019638530,0.0991563799137895,0.978925453896319,5.71604213160902,0.00584444499732144,4.34622614115380,205.842718056766,0.00232279820239211,0.714818981051436,14.3572871184989,0.00207116255606147,-9.51191877491359,0.961985951716582,0.492525625561131];
%             clear('initialPop')
            numberOfParametersToCalibrate = numel(lb);
            
            A = [zeros(1,16), 12, 1, 0, 0];
            b = 0;

        end
        
    end
end