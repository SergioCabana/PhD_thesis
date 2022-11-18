###################### PRIMERAS SIMULACIONES CON RASPASS (POSIBLEMENTE MAL, VE A LO SIGUIENTE) ####################

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

Hice ademas simulaciones a cenit 90deg, azimut 90deg cambiando el campo magnetico del polo sur y uno horizontal de 50 uT, con photonuclear off

Tengo ademas un set de 5 cascadas promediadas a distintas RASPASSHeights con interacciones fotonucleares y campo magnetico encendido o apagado,
para estudiar muones y sus fluctuaciones

Por ultimo tengo un conjunto de simulaciones promediando 5 cascadas a cenit 93,25 deg, azimut 90 deg ambiando RASPASSHeight y cambiando campos magneticos


Directorio en los nodos (input files, outputs):
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_proton    # primer conjunto
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_electron  # segundo conjunto
/home2/sergio.cabana/Test_RASPASS_ARENA22/Cascadas_foton     # tercer conjunto

######################### NUEVAS SIMULACIONES CON RASPASS ############################################
Creo que en todo lo anterior estaba exportando mal las tablas de datos, hay que indicar un Opt p para colocar planos perpendiculares al desarrollo
Aqui voy a repetir en cierto modo las simulaciones anteriores. Intentare ampliarlas de pase

# Primer conjunto

Foton de 1EeV, fijando RASPASSDistance para inyectar a 100km de altura, RHeight 15-20-25-30-35-40 km, RTimeShift 00ns

Todas las cascadas se hacen a cenit 93,25 deg para tener materia suficiente.

- Azimut 0deg, campo del polo sur: Promedio 5 cascadas con BothOn, GMOff, PNOff, BothOff
- Azimut 90deg, campo del polo sur: Promedio 5 cascadas con PNOff
- Azimut 270deg, campo del polo sur: promedio 5 cascadas con PNOff

- Azimut 0deg, campo horizontal: Promedio 5 cascadas con PNOff
- Azimut 90deg, campo horizontal: Promedio 5 cascadas con PNOff
- Azimut 270deg, campo horizontal: Promedio 5 cascadas con PNOff

Exporto tablas de e+-, e+, e- y si hay photonuclear, de mu+-, mu+, mu-

--- Aparentemente hay algun problema con azimut mayor a 90deg. El cambio en Opt p no ha sido muy importante (18/11/2022)

# Segundo conjunto:

Foton de 1EeV, fijando RASPASSDistance para inyectar en 100km de altura. RHeight 35, 36.5, 38, 40, 42 km (centrado
en valores grandes), RTimeShift 0ns

Todas las cascadas se hacen a cenit 93.25 deg

-Azimut 0deg, campo vertical y horizontal de 50uT
-Azimut 45deg, campo vertical y horizontal de 50uT
-Azimut 90deg, campo vertical y horizontal de 50 uT

Todo promediando 5 cascadas con PNOff, exporto desarrollo long. normal y unweighted de e+-, e+, e-


