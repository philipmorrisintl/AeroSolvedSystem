package ModelCases2p0
  model Human_Asgharian2001_Vt625ml "scenario: 
  - deposition in alveola (deep volume only)
  - alveola inlet volume = 10% of alveola total volume"
    extends
      Lung.WholeLungSinglePath.Simulators.BasicScenariosDeposition.Human_Asgharian2001(
      redeclare model WholeLungSP =
          Lung.WholeLungSinglePath.Models.SymmetricPathAirwaylDepositionNullInletVolAlv,
      system(redeclare package Medium =
            AeroSolvedSystem.Media.Mixtures.Water_air(CS=true,Nddc=22,Tm=Tm,pm=pm,Ym=Ym,Zm=Zm,Nm=Nm),
             useExtVolAlv = true,
             Nstart = Nm,
             Ystart = {0.0, 1-Nddc*zmm/1e5},
             Zstart = Zm/1e5),
      zmm=0.01/Nddc,
      Zm = [zmm,zmm,zmm,
            zmm,zmm,zmm,zmm,zmm,zmm,
            zmm,zmm,zmm,zmm,zmm,zmm,
            zmm,zmm,zmm,zmm,zmm,zmm,
            zmm;
            0,0,0,
            0,0,0,0,0,0,
            0,0,0,0,0,0,
            0,0,0,0,0,0,
            0] "Meadi Mixture substances liquid mass fractions for each diameter class",
      dmm = {5e-9,7e-9,9e-9,
             1e-8,1.3e-8,2e-8,3e-8,5e-8,7e-8,
             1e-7,1.3e-7,2e-7,3e-7,5e-7,7e-7,
             1e-6,1.3e-6,2e-6,3e-6,5e-6,7e-6,
             1e-5} "particle diameter for each diameter class",
      Vt = 0.625e-3,
      StartGenAl = 18 "start alveola region");
    extends AeroSolvedSystem.Icons.Validation;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=5,
        Interval=1.25,
        Tolerance=1e-06,
        __Dymola_Algorithm="Dassl"));
  end Human_Asgharian2001_Vt625ml;

  model Human_4regionsHeadML_Bf12_Vt625ml_Cycle_5_diam "scenario: 
  - deposition in alveola (deep volume only)
  - alveola inlet volume = 10% of alveola total volume"
    extends
      Lung.WholeLungSinglePath.Simulators.BasicScenariosDeposition.Human_4regionsHeadML(
      system(redeclare package Medium =
            AeroSolvedSystem.Media.Mixtures.Water_air(CS=true,Nddc=5,Tm=Tm,pm=pm,Ym=Ym,Zm=Zm,Nm=Nm),
             useExtVolAlv = true,
             Nstart = Nm,
             Ystart = {0.0, 1-Nddc*zmm/1e5},
             Zstart = Zm/1e5),
      zmm=0.01/Nddc,
      Zm = [zmm,zmm,zmm,zmm,zmm;
            0,0,0,0,0] "Meadi Mixture substances liquid mass fractions for each diameter class",
      dmm = {5e-9,
             1e-8,
             1e-7,
             1e-6,
             1e-5} "particle diameter for each diameter class",
      Bf = 12,
      Vt = 0.625e-3,
      StartGenAl = 18 "start alveola region",
      maskLungOut = false);
    extends AeroSolvedSystem.Icons.Validation;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=100,
        Interval=0.1,
        Tolerance=1e-04,
        __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput(events=false));
  end Human_4regionsHeadML_Bf12_Vt625ml_Cycle_5_diam;

  model Human_4regionsHeadML_Puff_Vp50ml_Vt625ml_2airCycle_5_diam "scenario: 
  - deposition in alveola (deep volume only)
  - alveola inlet volume = 10% of alveola total volume"
    extends
      Lung.WholeLungSinglePath.Simulators.BasicScenariosDeposition.Human_4regionsHeadML_Puff(
      system(redeclare package Medium =
            AeroSolvedSystem.Media.Mixtures.Water_air(CS=true,Nddc=5,Tm=Tm,pm=pm,Ym=Ym,Zm=Zm,Nm=Nm),
             useExtVolAlv = true,
             Nstart = Nm,
             Ystart = {0.0, 1-Nddc*zmm/1e5},
             Zstart = Zm/1e5),
      zmm=0.01/Nddc,
      Zm = [zmm,zmm,zmm,zmm,zmm;
            0,0,0,0,0] "Meadi Mixture substances liquid mass fractions for each diameter class",
      dmm = {5e-9,
             1e-8,
             1e-7,
             1e-6,
             1e-5} "particle diameter for each diameter class",
      Vt = 0.625e-3,
      Vp = 0.05e-3,
      StartGenAl = 18 "start alveola region",
      maskLungOut = false,
      puffingSinHoldBreakPrescP(AvLeakAir=1e-3));
    extends AeroSolvedSystem.Icons.Validation;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=60,
        Interval=0.1,
        Tolerance=1e-07,
        __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput(events=false));
  end Human_4regionsHeadML_Puff_Vp50ml_Vt625ml_2airCycle_5_diam;
  annotation (uses(Lung(version="4.0.0"), AeroSolvedSystem(version="4.0.0")));
end ModelCases2p0;
