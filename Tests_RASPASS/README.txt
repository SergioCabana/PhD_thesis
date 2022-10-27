###################### SIMULACIONES CON RASPASS ####################

# Primer conjunto

Condiciones: Proton de 1EeV, RDistance 1500km, RHeight 25-30-35-40 km, RTimeShift 0ns, cenit 93.25deg, azimut 0deg
Objetivo: reproducir resultados presentados por MAtias Tueros en ARENA 2022

Calculo en dominio temporal: Campo electrico observado en (0,0,RASPASSHeight)
Calculo en frecuencias: Cruz de antenas (plano yz) centrada en (0,0,RASPASSHeight), 100 antenas desde -7000m hasta 7000m del centro en cada brazo. Frecs [MHz]: 50, 100, 300, 500, 700, 1000
Altitud del suelo fijada a 0m (Site00), 350 Observing Levels, Saco desarrollo longitudinal de numero y energia de electrones y muones


# Segundo conjunto

Electron de 1EeV, RDistance fijada para inyeccion a 100km, RHeight 15-20-25-30-35-40 km, RTimeShift 0ns, cenit 93.25, azmiut 0
Simulaciones de desarrollo longitudinal, polo Sur

# Tercer conjunto

Foton de 1EeV, RDistance fijada para inyeccion a 100km, RHeight 15-20-25-30-35-40 km, RTimeShift 0ns, cenit 93.25, azmiut 0
Simulaciones de desarrollo longitudinal, polo Sur

Se hacen simulaciones con campo geomagnetico y procesos fotonucleares "activos" (BothOn), 
con el campo geomagnetico desactivado (GMOff), con los procesos fotonucleares descativados (PNOff),
y con ambos desactivados (BothOff) para estudiar efectos y diferencias

Directorio en los nodos (input files, outputs):
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_proton    # primer conjunto
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_electron  # segundo conjunto
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_foton     # tercer conjunto
