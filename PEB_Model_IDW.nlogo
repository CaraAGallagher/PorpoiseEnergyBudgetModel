;  PorpoiseModel version 2.0

; This is the source code for the individual-based model used for evaluating the
; effects of anthropogenic disturbances on the energetics and population dynamics of harbor porpoises
; in the inner Danish waters. The model builds on previous models of the effects of anthropogenic
; disturbances on porpoises in the IDW and North Sea. Please refer to the scientific publications for
; detailed documentation:

;      Gallagher, C.A., Grimm, V., Kyhn, L.A., Kinze, C., and Nabe-Nielsen, J. (2020).
;      "Movement and seasonal energetics mediate vulnerability to disturbance in marine
;      mammal populations". Submitted.

;      Nabe‐Nielsen, J., van Beest, F.M., Grimm, V., Sibly, R.M., Teilmann, J. and Thompson, P.M.,
;      (2018). "Predicting the impacts of anthropogenic disturbances on marine populations".
;      Conservation Letters, 11(5), p.e12563.

;      Nabe-Nielsen, J., Sibly, R.M., Tougaard, J., Teilmann, J. & Sveegaard, S. (2014)
;      "Effects of noise and by-catch on a Danish harbour porpoise population". Ecological
;      Modelling, 272, 242–251.

;      Nabe‐Nielsen, J., Tougaard, J., Teilmann, J., Lucke, K. and Forchhammer, M.C., (2013).
;      "How a simple adaptive foraging strategy can lead to emergent home ranges and
;      increased food intake". Oikos, 122(9), pp.1307-1316.

; The model was created as part of the projects:
; EFFECTS OF WIND FARMS ON HARBOUR PORPOISE BEHAVIOUR AND POPULATION DYNAMICS
; funded by Miljogruppen: By- og Landskabsstyrelsen, Energistyrelsen,
; DONG and Vattenfall A/S
; and
; BRIDGES AS BARRIERS PORPOISE MODELLING: DOES THE GREAT BELT BRIDGE HINDER
; MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT
; funded by Fehmern Belt A/S
; and
; The DEPONS project (www.depons.au.dk) funded by the offshore wind developers:
; Vattenfall, Forewind, ENECO Luchterduinen, Ørsted, and Scottish Power Renewables

; Copyright (C) 2016 (version 1.0), Jacob Nabe-Nielsen <jnn@bios.au.dk>
; Copyright (C) 2018 (version 1.1), Jacob Nabe-Nielsen <jnn@bios.au.dk>
; Copyright (C) 2020 (version 2.0), Jacob Nabe-Nielsen <jnn@bios.au.dk> & Cara Gallagher <cgallagher@bios.au.dk>

; This program is free software; you can redistribute it and/or modify it
; under the terms of the GNU General Public License version 2 and only
; version 2 as published by the Free Software Foundation.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

; The model was developed and tested using NetLogo version 6.1.0. Development ended 2020-03-04

; debug levels:
;   0   no debugging
;   1   debug porp-avoid-land / water depth problems
;   2   write turning angles before and after avoiding land
;   3   debugging turning angles -- esp. loop control
;   4   debugging attraction vector caused by reference memory
;   5   debugging porp-get-exp-food-val -- expected value of future food
;   6   debugging average energy level in 40x40 km blocks
;   7   debugging porpoise dispersal
;   8   debugging deterrence from ships and wind turbines
;   9   debugging reproduction
;   10  debugging energy budget procedures and energy related mortality
;   11  debugging energy intake procedures
;   11  debugging storage procedures




extensions [ gis profiler palette csv matrix] ;;; ADDED BY CARA: Palette csv Matrix
globals [
  time-step
  sim-day
  month                         ; using 30-day months
  quarter                       ; season of the year (used for updating food growth)
  year
  prev-month
  prev-quarter
  prev-year
  path                          ; path to directory for input/output, one dir for each area
  outfile                       ; simulation output
  bathy-data
  block-data                    ; block numbers 1-60 (40 x 40 km blocks, numbered by row, starting in upper left corner)
  block-centres-x               ; List of x coordinates for block centres
  block-centres-y               ; List of y coordinates for block centres
  block-values                  ; Expected relative quality of blocks, based on MEAN of Maxent values per block. Used for selecting dispersal direction
  list-of-dead-age              ; Age of death for all animals that die. Reset every year
  list-of-dead-day              ; Day of death for all animals that die. Reset every year
  disttocoast-data
  food-prob-data
  maxent-level-data
  mean-maxent-in-quarters       ; Standardized average maxent level in each quarter
  ; daily-surv-prob             ; Average daily survival in each porpoise age class. From Caswell et al 1998
  food-upd-interval             ; Number of days between updating amount of food. Updated daily in model >=1
  ;maxU                         ; Maximum utility of a patch, set to 1 here
  ;wind-turb-list               ; List of wind turbines. Each entry contains a list with id (char), x (num), y (num), and intensity (num)
  ;aarhus-odden-route           ; Positions on the ship route Aarhus-Odden. List of lists w two elements.
  ;great-belt-route             ; Positions on the ship T-route through Kattegat. List of lists w two elements.
  ;kattegat-sound-route         ; Positions on the ship T-route through the Sound. List of lists w two elements.
  age-distrib                   ; Number of porps in each age class (ages >=25 y aggregated)
  ;ship-ao-list                 ; List of ships in the Aarhus-Odden route
  ;ship-gb-list                 ; List of ships in the Great-Belt route
  ;ship-kat-list                ; List of ships in the Kattegat-Sound route
  ;turb-who-list                ;
  xllcorner
  yllcorner
  corr-logmov                   ; correlation in movement distance in CRW +
  corr-angle                    ; correlation in direction in CRW +
  deploy-x                      ; x-coord of porp 0
  deploy-y                      ; y-coord of porp 0
  turn-right                    ; random variable; turn right if = 1, left if = -1
  min-depth
  movlgt-list                   ; list of 100 move lengths for porpoise 0 (for histogram)
  angle-list                    ; list of 100 turning angles for porpoise 0 (for histogram)
  max-movlgt                    ; monotonously increasing, move length/100m
  mean-movlgt
  memory-max                    ; Maximum number of half-hour steps the amount of food can be remembered (120 steps is 2.5 days)
  inertia-const                 ; A, set to 0.001 like in ms
  ref-mem-strength-list-fixed   ; replaces ref-mem-strength-list, uses rR = 0.10
  work-mem-strength-list-fixed  ; replaces work-mem-strength-list, uses rW = 0.20
  use-exp-food-val?              ; get more attracted to the CRW path if food was found recently
  vt                            ; resultant attraction vector, resulting from reference memory of food availability (model >=2)
  CRW-contrib                   ; length of vector pointing in direction predicted by CRW
  ; MR-contrib                  ; length of vector pointing in direction of remembered food

  ;;; ********************* ADDED BY CARA ********************** ;;;

  ;;; Slider params from Nabe-Nielsen et al. 2014 not changed in this version of the model:--;;;
  ; dispersion params
  min-disp-depth          ; Minimum dispersal depth in m
  n-disp-targets          ; Number of dispersal targets
  mean-disp-dist          ; Mean dispersal distance in km / 30 min
  min-dist-to-target      ; Minimum distance to target in km

  ; food growth
  maxU                    ; Maximum utility of a patch, set to 1 here
  food-growth-rate        ; Food patch regrowth rate (rU)
  ;;;-----------------------------------------------------------------------------------------;;;

  ; porpoise globals
  list-of-dead-age-calves ; Age of death for all calves that die. Reset every year
  list-of-dead-day-calves ; Day of death for all calves that die. Reset every year
  IR-temp-mod             ; seasonal modifier of intake rate and mean storage level based on water temperature
  IR-coef                 ; maximum ingestion rate scaling coefficient
  satiation-c             ; satiation constant, unitless
  AE-food                 ; assimilation efficiency of food, %
  IR-to-EA                ; conversion from ingested food to energy in joules
  prop-eff                ; propulsive efficiency of swimming, %
  age-of-maturity         ; age when becoming receptive, in years
  pregnancy-rate          ; pregnancy rate of the population
  repro-min-SL            ; Reproductive storage level threshold, %
  calf-idl-SL             ; Calf ideal storage level, %
  max-mass-f              ; neonate birth size, kg
  f-growth-c              ; fetal growth constant, unitless
  percent-lip-f           ; percent lipid fetus, %
  percent-pro-f           ; percent protein of the fetus, %
  lact-eff                ; lactation efficiency, %
  t-gest                  ; gestation time, days
  t-nurs                  ; nursing time, days
  ED-lip                  ; energy density of lipid, J kg-1
  ED-pro                  ; energy density of protein, J kg-1
  dens-blub               ; density of lipid, kg cm-3
  x-surv-prob-const       ; survival probability constant

  ; environmental globals
  temp-data               ; temperature maps
  salinity-data           ; salinity maps
  density-matrix          ; lookup table of density values for given temperature and salinity values
  row-dens
  col-dens
  conduct-matrix          ; lookup table of water conductivity values for given temperature and salinity values
  row-conduct
  col-conduct
  ceote-matrix            ; lookup table of coefficient of thermal expansion values for given temperature and salinity values
  row-ceote
  col-ceote
  sp-heat-matrix          ; lookup table of specific heat values for given temperature and salinity values
  row-sp-heat
  col-sp-heat
  dynamic-vis-matrix      ; lookup table of water dynamic viscocity values for given temperature and salinity values
  row-dyn-vis
  col-dyn-vis
  gravity                 ; gravitational constant in m s-1
  mean-temp               ; average monthly temperature
  water-patches           ; patches of water agentset


  ; added for scenarios
  no-surveys
  ss-ship-order
  n-calf-lost
  daily-pop-list
  fourty-km-list
  storage-level-deter-list
  mean-calf-mass-deter-list
  mean-juv-mass-deter-list
  mean-calf-mass-list
  mean-juv-mass-list
  fourty-km-out
  storage-level-out
  mean-calf-mass-out
  mean-juv-mass-out
  tot-calf-loss-out
  ]


breed [ porps porp ]
breed [ ss-ships ss-ship ]   ; ADDED BY CARA
;breed [ turbs turb ]
;breed [ ao-ships ao-ship ]
;breed [ gb-ships gb-ship ]
;breed [ kat-ships kat-ship ]


patches-own [
  bathymetry
  block
  disttocoast            ; Distance to coast (in m)
  food-prob              ; randomly distributed patches, 1 cell large (prob = 1 inside patches, 0 outside)
  maxent-level           ; growth rate for the food in the patches (from quarterly MAXENT predictions)
  food-level             ; current amount of food in cell (zero outside food patches)

  ;;; ********************* ADDED BY CARA ********************** ;;;
  temp-w                       ; temperature of the water in C
  salinity-w                   ; water salinity in g kg-1
  density-w                    ; water density in kg m-3
  conductivity-w               ; water conductivity in W C-1 m-1
  coef-thermal-expansion-w     ; coefficient of thermal expansion in C-1
  Pr-w                         ; prandtl number, unitless
  kin-visc-w                   ; kinematic viscocity of water in m2 s-1
  deep-water?                   ; true if distance from coast greater than max movement length
  ]
porps-own [
  lgth                       ; Length in cm
  weight                     ; Weight in kg
  age                        ; Age in years (remember, 360 days per year)
 ; age-of-maturity           ; Age when becoming receptive, in years ; CHANGED BY CARA
  pregnancy-status           ; 0 (unable to mate, young/low energy); 1 (unable to mate, pregnant); 2 (ready to mate)
  mating-day                 ; Day of year. Most mating occurs in August (Lockyer 2003)
  ds-mating                  ; Days since mating. -99 if not pregnant
  dsg-birth                  ; Days since giving birth. -99 if not with lactating calf
  with-lact-calf?            ; true/false, with lactating calf
  ;energy-level              ; Porpoises get energy by eating and loose energy by moving.       ;REMOVED BY CARA
  ;energy-level-sum          ; Sum of energy levels. Reset to 0 every day                       ;REMOVED BY CARA
  ;energy-level-daily        ; List with average energy for last ten days. Latest days first.   ;REMOVED BY CARA
  disp-type                  ; Disperse away from low-energy area. 0 if not dispersing, 1 if dispersing far, 2 if...
  disp-target                ; List with x and y coord of the patch that the porp attempts to disperse to (not UTM)
  prev-angle                 ; Last turning angle (not the heading!)
  pres-angle                 ; Present turning angle
  prev-logmov                ; Previous Log10 (move length [measured in 100-m steps])
  pres-logmov                ; Present Log10 (move length [measured in 100-m steps])
  enough-water-ahead?         ; Turn to avoid land if false
  pos-list                   ; Coordinates of previous positions -- one per 30 min, latest positions in left end (= item 0)
  pos-list-daily             ; Coordinates of previous 10 daily positions -- daily, corresponding to energy-level-daily
  ; Vars added in model 2
  ; ref-mem-strength-list    ; Memory decay with time (logistic decrease, function of rR)
  ; work-mem-strength-list   ; Memory of poor quality patches, decays with time (logistic decrease, function of rW)
  ; work-mem-updated         ; Whether a list of working memory has been updated
  stored-util-list           ; Up-t ; Remembered feeding success (after memory decay) -- latest positions left
  deter-vt                   ; Vector (= list) determining which direction a porp is deterled from wind turbines and ships, and how much
  deter-strength
  deter-time-left
  VE-total                   ; total value of food expected to be found in the future


  ;;; ********************* ADDED BY CARA ********************** ;;;
  storage-level           ; Porpoises get energy by eating and loose energy by moving.
  storage-level-sum       ; Sum of energy levels. Reset to 0 every day
  storage-level-daily     ; List with average energy for last ten days. Latest days first.


  age-days                          ; age in days
  swim-speed                        ; swimming speed in m s-1
  e-storage                         ; mobilizable energy available in storage in J
  e-repo-min                        ; minimum storage energy needed for reproduction in J
  mass-struct                       ; structural mass (non-blubber) in kg
  mass-struct-calf                  ; structural mass of the calf in kg
  v-blub                            ; percent body mass represented by blubber
  v-blub-mean                       ; mean value of blubber, representing porpoises in good body condition in cm3
  v-blub-min                        ; minimum value of blubber before starving in cm3
  v-blub-repro                      ; minimum blubber volume needed for reproduction in cm3
  v-blub-calf                       ; the calf's volume of blubber in cm3
  v-blub-calf-idl                   ; the ideal volume of blubber for the calf in cm3
  SL-mean                           ; mean storage level as percentage of body mass
  sex-calf                          ; sex of calf
  mass-f                            ; mass of the fetus in kg
  mass-calf                            ; mass of the calf in kg
  m-str-k                           ; mass growth constant k, unitless
  m-str-inf                         ; asymptotic mass in kg
  lgth-calf                         ; length of the calf in m
  lgth-inf                          ; asymptotic length in cm
  lgth-0                            ; minimum length  in cm
  lgth-k                            ; length estimator k, unitless
  surface-area                      ; surface area of porpoise in m2

; Energy intake
  e-assim                           ; assimliated energy in joules
  e-assim-grow                      ; e-assim available for growth (split between growth and blubber)
  e-assim-grow-calf                 ; e-assim available for calf growth
  e-assim-calf                      ; energy assimilated by feeding calves
  IR-record                         ; maximum ingestion rate in kg
  IR-record-calf                    ; ingestion rate of the dependant calf in kg

; Maintenance
  B0                                ; basal metabolic rate normalization constant, unitless

; Thermoregulation
  for-convec-scal-coef-list         ; forced convection scaling coefficient in W m−2 °C−1
  hC-low-lim-coef-list              ; heat transfer coefficient lower limit W m−2 °C−1
  T-core                            ; core temperature of the animal in °C
  kB                                ; thermal conductivity of blood free blubber in W m−1 °C−1

; Locomotion
  lambda                            ; ratio of active to passive drag

; Growth and Reproduction
  max-grow                          ; maximum growth as determined by gompertz relationship in kg timestep-1
  max-grow-calf                     ; maximum growth of the calf in kg timestep-1
  struct-mass-perc-pro              ; structural mass percent protein
  struct-mass-perc-lip              ; structural mass percent lipid
  DE-lip                            ; deposition efficiency of lipid as a %
  DE-pro                            ; deposition efficiency of protein as a %
  ED-lean-mass                      ; energy density of lean mass in J kg-1
  DE-lean-mass                      ; deposition efficiency of lean mass as a %
  wean-scale-fact                   ; weaning scale factor as a %
  preg-chance

; Storage
  perc-lip-blub                     ; percent lipid of blubber
  ; Storage: blubber depth sites
  axillary-D                        ; axillary dorsal site blubber depth in cm
  axillary-L                        ; axillary lateral site blubber depth in cm
  axillary-V                        ; axillary ventral site blubber depth in cm
  CrIDF-D                           ; cranial insertion of the dorsal fin dorsal site blubber depth in cm
  CrIDF-L                           ; cranial insertion of the lateral fin dorsal site blubber depth in cm
  CrIDF-V                           ; cranial insertion of the ventral fin dorsal site blubber depth in cm
  CaIDF-D                           ; caudal insertion of the dorsal fin dorsal site blubber depth in cm
  CaIDF-L                           ; caudal insertion of the lateral fin dorsal site blubber depth in cm
  CaIDF-V                           ; caudal insertion of the ventral fin dorsal site blubber depth in cm

  ; Storage: cone attributes used for thermo submodel
  c2                                ; list containing average height, length, and blubber depth for cone 2
  c3                                ; list containing average height, length, and blubber depth for cone 3

; Energy budget outputs
  m-BMR                             ; basal metabolic rate in J timestep-1
  m-BMR-calf                        ; energy needed for calf BMR in J timestep-1
  m-thermo                          ; metabolic cost of thermoregulation in J timestep-1
  m-thermo-calf                     ; energy needed to for calf thermoregulation in J timestep-1
  m-loco                            ; metabolic cost of locomotion in J timestep-1
  m-loco-ineff-pts                  ; locomotion inefficiency from the previous timestep in J timestep-1
  m-growth                          ; metabolic cost of growth in J timestep-1
  m-growth-f                        ; metabolic cost of growth for the fetus in J timestep-1
  m-growth-calf                     ; energy needed for calf growth in J timestep-1
  e-heat-gest                       ; metabolic cost of the heat of gestation in J timestep-1
  m-preg                            ; metabolic cost of pregnancy in J timestep-1
  m-blub-calf                       ; energy needed to maintain calf's blubber stores in J timestep-1
  blub-calf                         ; amount of blubber to be allocated to calves in cm3
  e-calf                            ; total calf costs for lactating calves in J timestep-1
  m-lact                            ; metabolic cost of lactation in J timestep-1
  m-lact-real                       ; cost actually spent after allocation procedures in J timestep-1
  m-tot                             ; total metabolism in J timestep-1
  vital-costs                       ; vital costs of calf to be covered by mother in J
  vital-costs-calf                  ; vital costs of calf to be covered by calf in J

; counters
  abortion-count                    ; abortion counter

  ; For testing:
  daily-food
  food-intake-list

 ]

ss-ships-own [ ; ADDED BY CARA
  seismic-survey-route   ; Positions on the seismic survey route. List of lists w two elements.
  ss-order
  start-ts
  end-ts
  route
  id
  ship-speed
  impact
  past-loc
  next-loc
  TL-B
  TL-a
]



to landsc-setup
  ; Note that setting the coordinate system here is optional, as
  ; long as all of your datasets use the same coordinate system.
  ; gis:load-coordinate-system (word "data/" projection ".prj")
  ; Load datasets:
  set path word "environmentals/raster-data/" area      ; note that "/" may be a problem on Windows
  set path word path "/"
  set bathy-data gis:load-dataset word path "bathy.asc"
  set food-prob-data gis:load-dataset word path "patches.asc"
  set maxent-level-data gis:load-dataset word path "quarter1.asc"
  set block-data gis:load-dataset word path "blocks.asc"
  set block-centres-x ( list 50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550 )
  set block-centres-y ( list 950 950 950 950 950 950 850 850 850 850 850 850 750 750 750 750 750 750 650 650 650 650 650 650 550 550 550 550 550 550 450 450 450 450 450 450 350 350 350 350 350 350 250 250 250 250 250 250 150 150 150 150 150 150  50  50  50  50  50  50 )
  set disttocoast-data gis:load-dataset word path "disttocoast.asc"

  ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;

  set temp-data gis:load-dataset word path "TemperatureData/TempData_Month1.asc"
  set salinity-data gis:load-dataset word path "SalinityData/SalData_Month1.asc"

  ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

  ; Init landscape-related variables
  set food-upd-interval 1
  ;set maxU 1
  set time-step 0
  set sim-day 0
  set year 0
  set month 1
  set quarter 1

  ; Use Kattegat-area throughout
  set xllcorner 529473
  set yllcorner 5972242

  ; Set the world envelope to the union of all of our dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of
    (gis:envelope-of bathy-data)
    (gis:envelope-of food-prob-data)
    (gis:envelope-of maxent-level-data)

    ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;

    (gis:envelope-of temp-data)
    (gis:envelope-of salinity-data)

   ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;
    )

  ; This is the preferred way of copying values from a raster dataset
  ; into a patch variable: in one step, using gis:apply-raster.
  gis:apply-raster bathy-data bathymetry
  gis:apply-raster food-prob-data food-prob
  gis:apply-raster maxent-level-data maxent-level
  gis:apply-raster block-data block
  gis:apply-raster disttocoast-data disttocoast

  ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;

  gis:apply-raster temp-data temp-w
  gis:apply-raster salinity-data salinity-w

  ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

  ; set amount of food -- if there is a chance that food is present
  set mean-maxent-in-quarters (list 0.515686364223653 0.888541219760357 0.841346010536882 1) ; setting global variable here
  ask patches [ ifelse not ( maxent-level >= 0 ) [ set maxent-level -1 ] [ set maxent-level maxent-level ] ]  ; ### new check in PorpoiseModel version 1.1
  ask patches [ ifelse food-prob > 0 [ set food-level maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters ] [ set food-level 0 ] ]
  ;ask patches [ ifelse food-prob > 0 [ set food-level maxU * 1 / item (quarter - 1) mean-maxent-in-quarters ] [ set food-level 0 ] ]
  ask patches [ if food-prob > 0 and food-level < 0.01 [ set food-level 0.01 ]]

 ; ask patches [
 ;   ifelse food-prob > 0 and maxent-level > 0 [ set food-level (maxU *  (sqrt maxent-level / 0.92))  / item (quarter - 1) mean-maxent-in-quarters ] [ set food-level 0 ]  ; CHANGED BY CARA; as in PorpoiseModel version 1.1, but changed to include maxent level from beginning to ensure relative rates food levels are correct
 ;   if food-prob > 0 and food-level < 0.01 [ set food-level 0.01 ]
 ; ]


  ; Next, handle NaNs in distance to coast data. Test introduced in PorpoiseModel version 1.1
  ask patches [ ifelse not ( disttocoast >= 0 ) [ set disttocoast 0 ] [ set disttocoast disttocoast ] ]
  ask patches [ ifelse not ( bathymetry >= 0 ) [ set bathymetry 0 ] [ set bathymetry bathymetry ] ]
  ask patches [ ifelse not ( block >= 0 ) [ set block -1 ] [ set block block ] ]
  ask patches [ ifelse not ( food-prob >= 0 ) [ set food-prob 0 ] [ set food-prob food-prob ] ]

  ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;

  ; set patch temperature and salinity values
  ask patches [ ifelse not ( temp-w >= 0 ) [ set temp-w -1 ] [set temp-w temp-w ] ]
  ask patches [ ifelse not ( salinity-w >= 0 ) [ set salinity-w -1 ] [ set salinity-w salinity-w ] ]
  set water-patches patches with [temp-w > 0 and salinity-w > 0 and bathymetry > 0]

  ; if porpoises cannot possibly reach the coast in one timestep from that location then set it to be deep water
  ask patches [ ifelse ((disttocoast / 400) >= 15.25) [ set deep-water? true ] [ set deep-water? false ]]

  ; import thermophysical properties of seawater lookup tables
  ; density
  let density csv:from-file "environmentals/thermophysical properties of SW tables/DensityTable.csv"
  set density-matrix matrix:from-row-list density
  set row-dens matrix:get-row density-matrix 0
  set col-dens matrix:get-column density-matrix 0

  ; thermal conductivity of water
  let conduct csv:from-file "environmentals/thermophysical properties of SW tables/ThermalConductivity.csv"
  set conduct-matrix matrix:from-row-list conduct
  set row-conduct matrix:get-row conduct-matrix 0
  set col-conduct matrix:get-column conduct-matrix 0

  ; coef of thermal expansion
  let ceote csv:from-file "environmentals/thermophysical properties of SW tables/CoefOfThermalExpansionTable.csv"
  set ceote-matrix matrix:from-row-list ceote
  set row-ceote matrix:get-row ceote-matrix 0
  set col-ceote matrix:get-column ceote-matrix 0

  ; specific heat of water
  let sp-heat csv:from-file "environmentals/thermophysical properties of SW tables/SpecificHeat.csv"
  set sp-heat-matrix matrix:from-row-list sp-heat
  set row-sp-heat matrix:get-row sp-heat-matrix 0
  set col-sp-heat matrix:get-column sp-heat-matrix 0

  ; dynamic viscocity
  let dyn-viscosity csv:from-file "environmentals/thermophysical properties of SW tables/DynamicVis.csv"
  set dynamic-vis-matrix matrix:from-row-list dyn-viscosity
  set row-dyn-vis matrix:get-row dynamic-vis-matrix 0
  set col-dyn-vis matrix:get-column dynamic-vis-matrix 0

  ; patches then pull their thermophysical properties based on their temperature and salinity values
  ask water-patches [
      thermo-prop-update
      ]

  ; estimate intake rate modifier based on mean water temp
  set mean-temp mean [temp-w] of water-patches
  let xxxxx ((1 / 6.3661) * mean-temp) * (180 / pi)
  let yyyyy (cos xxxxx)
  set IR-temp-mod ((1 / 5)* yyyyy + 1)

  ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

  landsc-display

  print ""
  let tmp-t word "Setup time: " timer
  print word tmp-t " sec"
end


; Turbine code removed for this application::: CHANGED BY CARA

to file-noise-debug ; setup file and write; Turbines in line have ycor = 706
    let outfile-for-debug "output/Homogeneous/file-noise-debug.txt" ; overwrite existing data
    if (not file-exists? outfile-for-debug) [
      file-open outfile-for-debug
      print ( "Writing to output/Homogeneous/file-noise-debug.txt" )
      file-print (" deterr-coef std-deter-dist id time-step utm-x utm-y") ; header line (space-separated).
      file-close
    ]
    file-open outfile-for-debug
    ask porps [
      ;if ( ycor < (706 + 5) and ycor > (706 - 5)) [
      if ( ycor < (706 + 20) and ycor > (706 - 20)) [
        file-write deterrence-coeff
        file-write std-deterrence-dist
        file-write who
        file-write time-step
        file-write xcor
        file-write ycor
        file-print ""    ; followed by CR
      ]
    ]
    file-close
end ; end file-noise-debug


 to ships-import  ; get ship impacts and speeds [km / h] in all routes. Separate lists for different routes: ; CHANGED BY CARA
  ;; aarhus-odden route: ship-ao-list read from file
;  if (incl-ships?) [
;    let ship-ao "ships/Aarhus-Odden.txt"
;    let ship-gb "ships/Great-Belt.txt"
;    let ship-kat "ships/Kattegat-Sound.txt"
;    let line-txt ""
;    let line-lst ( list " "  " "  " " )
;    set ship-ao-list ( list line-lst ); a list of line-lst elements
;    set ship-gb-list ( list line-lst ); a list of line-lst elements
;    set ship-kat-list ( list line-lst ); a list of line-lst elements
;    let line-lgt 0
;    if (not file-exists? ship-ao) [
;      print word "Cannot load file: " ship-ao
;      stop
;    ]
;    if (not file-exists? ship-gb) [
;      print word "Cannot load file: " ship-gb
;      stop
;    ]
;    if (not file-exists? ship-kat) [
;      print word "Cannot load file: " ship-kat
;      stop
;    ]
    ;; One route at a time -- close file afterwards:
;    file-open ship-ao
;    let i 0  ; chars
;    let j 0  ; word nbr
;    let k 0  ; line nbr
;    let one-char ""
;    let one-word ""
;    let in-word? true
;    while [ not file-at-end? ] [
;      set line-txt file-read-line
;      set line-txt word line-txt " "
;      set line-lgt length line-txt
;      set i 0
;      set j 0
;      set in-word? true
;      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
;        set one-char substring line-txt i (i + 1)
;        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word?) [  ; that is, prev. char was part of a word, current isn't
;          set in-word? false
;          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
;          set line-lst replace-item j line-lst one-word
;          set j j + 1
;        ]
;        if ((one-char = "\t" or one-char = " ") and (not in-word?)) [
;          set one-word ""
;        ]
;        if (one-char != "\t" and one-char != " ") [
;          set in-word? true
;          set one-word word one-word one-char
;        ]
;        set i i + 1
;      ]
;      if (k > 0) [ set ship-ao-list lput line-lst ship-ao-list ]
;      set k k + 1
;    ]
;    file-close
    ; next route:
 ;   file-open ship-gb
 ;   set i 0  ; chars
 ;   set j 0  ; word nbr
 ;   set k 0  ; line nbr
 ;   set one-char ""
 ;   set one-word ""
 ;   set in-word? true
 ;   while [ not file-at-end? ] [
 ;     set line-txt file-read-line
 ;     set line-txt word line-txt " "
 ;     set line-lgt length line-txt
 ;     set i 0
 ;     set j 0
 ;     set in-word? true
 ;     while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
 ;       set one-char substring line-txt i (i + 1)
 ;       if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word?) [  ; that is, prev. char was part of a word, current isn't
 ;         set in-word? false
 ;         if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
 ;         set line-lst replace-item j line-lst one-word
 ;         set j j + 1
 ;       ]
 ;       if ((one-char = "\t" or one-char = " ") and (not in-word?)) [
 ;         set one-word ""
 ;       ]
 ;       if (one-char != "\t" and one-char != " ") [
 ;         set in-word? true
 ;         set one-word word one-word one-char
 ;       ]
 ;       set i i + 1
 ;     ]
 ;     if (k > 0) [ set ship-gb-list lput line-lst ship-gb-list ]
 ;     set k k + 1
 ;   ]
 ;   file-close
   ;; next route:
 ;   file-open ship-kat
 ;   set i 0  ; chars
 ;   set j 0  ; word nbr
 ;   set k 0  ; line nbr
 ;   set one-char ""
 ;   set one-word ""
 ;   set in-word? true
 ;   while [ not file-at-end? ] [
;      set line-txt file-read-line
;      set line-txt word line-txt " "
;      set line-lgt length line-txt
;      set i 0
;      set j 0
;      set in-word? true
;      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
;        set one-char substring line-txt i (i + 1)
;        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word?) [  ; that is, prev. char was part of a word, current isn't
;          set in-word? false
;          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
;          set line-lst replace-item j line-lst one-word
;          set j j + 1
;        ]
;        if ((one-char = "\t" or one-char = " ") and (not in-word?)) [
;          set one-word ""
;        ]
;        if (one-char != "\t" and one-char != " ") [
;          set in-word? true
;          set one-word word one-word one-char
;        ]
;        set i i + 1
;      ]
;      if (k > 0) [ set ship-kat-list lput line-lst ship-kat-list ]
;      set k k + 1
;    ]
;    file-close

;  ]
 end ; end ships-import


to ships-setup ; CHANGED BY CARA
  ; setup seismic survey vessels
  let n-ss-ships no-surveys
  create-ss-ships n-ss-ships [
      ; assign tracks using a list of buoy coordinates:
      ; Omø location
      ifelse (item 0 ss-ship-order) = 1 [ set seismic-survey-route (list (list 297 343)	(list 290.5 320.5) (list 287 317.5)	(list 255 325) (list 248 330.5)	(list 288 321) (list 292 338)	(list 289 341) (list 283.5 317.5)	(list 291.5 323.5)	(list 240 336)	(list 237.5 340)	(list 294.5 326.5)	(list 295.5 330)	(list 238 345)	(list 238.5 348)	(list 298.5 333.5)	(list 299 336.5)	(list 239 352)	(list 243 369.5)	(list 264.5 364)	(list 264 360.5)	(list 242 365.5)	(list 241.5 362.5)	(list 263.5 357)	(list 263.5 353.5)	(list 240.5 359)	(list 240 355.5)	(list 285.5 343.5)	(list 280 319)	(list 277 319.5)	(list 282.5 346)	(list 279.5 348.5)	(list 273.5 323.5)	(list 271 324)	(list 276 348.5)	(list 272.5 350)	(list 265.5 319.5)	(list 262 320.5)	(list 268.5 349.5)	(list 266 350.5)	(list 259 323)	(list 256 324)	(list 267 373)	(list 263.5 373.5)	(list 252.5 326.5)	(list 249.5 329)	(list 260.5 378.5)	(list 257.5 380.5)	(list 246.5 331.5)	(list 243.5 333.5)	(list 254.5 382.5)	(list 250.5 382)	(list 257.5 380.5)	(list 261.5 376.5)	(list 248.5 379)	(list 247.5 375.5)	(list 266 371)	(list 265.5 367.5)	(list 245 372)	(list 237.5 340)	(list 240 336)	(list 250.5 385.5))
      set route "seismic-survey-Omø"]
      ; Sejerø Bugt location
      [ifelse (item 0 ss-ship-order) = 2 [set seismic-survey-route (list (list 312.5 572) (list 304.5 549) (list 300.5 547) (list 269.5 557) (list 263 563) (list 302.5 550) (list 307.5 568) (list 305 571) (list 297.5 547) (list 305.5 552) (list 255.5 570) (list 253.5 574) (list 309 555) (list 310.5 558) (list 254 578) (list 255 581) (list 313.5 561.5) (list 314.5 564) (list 256 585) (list 261.5 602) (list 282 595) (list 281.5 591.5) (list 260 598.5) (list 259.5 595) (list 280.5 588) (list 280.5 584.5) (list 258 591.5) (list 257 588.5) (list 301.5 573.5) (list 294 549) (list 291 549.5) (list 299 575.5) (list 296 578) (list 288 554) (list 285.5 554.5) (list 292.5 578) (list 289 580) (list 279.5 550.5) (list 276.5 551.5) (list 285.5 580) (list 282.5 581.5) (list 273.5 554.5) (list 270.5 555.5) (list 285.5 603.5) (list 282.5 604.5) (list 267.5 559) (list 264.5 561) (list 279.5 609.5) (list 276.5 612) (list 261.5 564) (list 258.5 567.5) (list 273.5 614) (list 270 614) (list 276.5 612) (list 280.5 607.5) (list 267.5 611) (list 266.5 608) (list 284.5 601.5) (list 283.5 598.5) (list 263.5 604.5) (list 253.5 574) (list 255.5 570) (list 270.5 618))
      set route "seismic-survey-Sejerø"]
      ; Samsø location
      [set seismic-survey-route (list (list 153.5 513)	(list 152 490)	(list 148.5 486.5)	(list 116 487)	(list 108 491)	(list 149.5 490)	(list 149.5 507.5)	(list 146.5 509.5)	(list 145.5 485.5)	(list 152 493)	(list 99 495)	(list 96 498.5)	(list 154.5 496.5)	(list 154.5 500.5)	(list 95 502)	(list 95 505.5)	(list 157 504)	(list 157 507)	(list 95 509.5)	(list 95.5 527.5)	(list 117.5 526.5)	(list 118 523)	(list 95.5 523.5)	(list 95.5 520)	(list 118 519)	(list 119 516)	(list 95.5 516.5)	(list 95 513)	(list 142 511)	(list 142 486.5)	(list 138.5 486.5)	(list 139 512.5)	(list 135.5 514)	(list 134.5 489.5)	(list 131.5 489.5)	(list 132 513.5)	(list 128.5 514)	(list 127.5 484)	(list 124 484)	(list 124.5 513.5)	(list 121.5 513.5)	(list 120.5 486)	(list 117 486.5)	(list 118.5 535.5)	(list 115 535.5)	(list 113 488.5)	(list 110 490)	(list 111.5 540)	(list 107.5 541)	(list 106.5 495)	(list 103 493.5)	(list 104 542.5)	(list 100.5 541.5)	(list 107.5 541)	(list 112.5 538)	(list 99 538)	(list 99 534.5)	(list 118 533.5)	(list 118 530)	(list 97 530.5)	(list 96 498.5)	(list 99 495)	(list 100.5 545))
       set route "seismic-survey-Samsø"]]

      set past-loc 0 ; start at pos 1 in survey
      set shape "boat top"
      set color black
      set size 25
      setxy (item 0 (item 0 seismic-survey-route)) (item 1 (item 0 seismic-survey-route)) ; set position at first coordinate
      facexy (item 0 (item next-loc seismic-survey-route)) (item 1 (item next-loc seismic-survey-route)) ] ; face second coordinate

  ask ss-ships [
      set id 1
      set ship-speed 0.873  ; speed allowing for track completion in 30 days
      set impact source-level ; source level in dBpeak-peak from Kyhn et al. 2019

      set start-ts (((sc-year - 1) * 360 * 48) + ((sc-month - 1) * (30 * 48))) ; starting timestep
      set end-ts (start-ts + (1 * 30 * 48)) ; ending timestep

      ; setup transmission loss parameters for each survey site:
      if route = "seismic-survey-Omø" [
        set TL-B 16.7
        set TL-a 0.00086
    ]
      if route = "seismic-survey-Sejerø" [
        set TL-B 18.7
        set TL-a 0.000343
    ]
       if route = "seismic-survey-Samsø" [
        set TL-B 17.7
        set TL-a 0.0006
    ]
  ]
end

to ships-move ; CHANGED BY CARA
  ; ships move along designated routes by following buoy coordinates
  let the-route "NA"
  let mov-lgt 0

    ask ss-ships [
    if time-step >= start-ts and time-step <= end-ts [
    pen-down
    set mov-lgt ship-speed
    set the-route seismic-survey-route
    ; approaching next location on route?
    if distancexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route)) < mov-lgt [
      ; check wheter position numbers should increase, or if ships should turn around:
      let tmp next-loc
      ifelse (next-loc > past-loc and next-loc < (length(the-route) - 1) ) [ set next-loc next-loc + 1 ] [
      set next-loc next-loc - 1 ]
      set past-loc tmp
      if  (next-loc < 0) [ set next-loc 1 ]
      facexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route))
    ]
    fd mov-lgt
  ]
    if scenario != 0 [if time-step > end-ts [ die ]]
      ]
end


to ships-deter-porps ; CHANGED BY CARA
  ; ships produce sound and scare porpoises
  let curr-deter 0 ; current amount of deterring
  let deter-dist std-deterrence-dist / 400 ; number of grid-cells where a ship will potentially affect a porpoise
  let ship-pos list -9 -9
  let porp-pos list -9 -9
  let dist-to-ship -9
  let rec-level 0
  ask ss-ships [
  if time-step >= start-ts and time-step <= end-ts [
  set ship-pos list ([xcor] of self) ([ycor] of self)
  ask porps in-radius deter-dist [  ; i.e. porps < deter-dist away
      set porp-pos list ([xcor] of self) ([ycor] of self)                  ; porpoise position
      set dist-to-ship distancexy (item 0 ship-pos) (item 1 ship-pos)      ; calculate distance from ship
      let r (dist-to-ship * 400)                                           ; meters away from ship

      ; check if land between ship and porp
      let head heading
      let dist 1
      let land-ahead? false
      facexy (item 0 ship-pos) (item 1 ship-pos)
      while [dist <= dist-to-ship] [
      let p patch-ahead dist
      if [disttocoast] of p = 0 [set land-ahead? true]
      set dist dist + 1
      ]
      set heading head

      ifelse land-ahead? = false
      [ set rec-level ([impact] of myself) - ([TL-B] of myself * (log r 10)) - ([TL-a] of myself * r) ] ; received sound level based on RL = SL - (Beta log10(r) + alpha(r))
      [ set rec-level 0 ]                                                                               ; if land barrier recieved level = 0

      ifelse rec-level > resp-level
      [ set curr-deter (rec-level - resp-level) ] ; deterring-strength decreases logistically with distance as in DEPONS - deterrence coefficient in porp-std-move
      [ set curr-deter 0 ]

      if (deter-strength < curr-deter ) [ ; become deterred if not already more scared
        set deter-strength curr-deter
        set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 ship-pos)))   ; vector pointing away from turbine
        set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 ship-pos)))
        set deter-time-left deter-time ; how long to remain affected
      ]
      if (debug = 8) [
       ; if land-ahead? = false [
        if deter-strength > 0 [
        set color 105
        ;type word "who: " who
        ;type word " dist-to-ship " myself
        ;type word ": " dist-to-ship
        ;print word " received-level: " rec-level
         ]
          ]

      ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
      set disp-type 0
     ]
  ]]

end ; end ships-deter-porps


to update-block-values
  ; Update relative values for blocks 1-60 every quarter, depending on location -- based on MEAN of maxent values for block
  let block-val-homo ( list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 )
  let block-val-kat-1 ( list 0 0.007 0.0042 0.0022 8e-04 0 0 0.0454 0.0522 0.0176 0.0185 0 0 0.0458 0.0348 0.0062 0.0182 0.4446 0 0.2898 0.1013 0.0769 0.126 0.4615 0.4439 0.5104 0.5741 0.4593 0.6447 0.6498 0.8685 0.5063 0.971 0.1697 0.27 0.1261 0.9438 0.7346 1 0.3994 0.1576 0.2355 0.7718 0.9082 0.8524 0.2161 0.4916 0.2975 0.8397 0.7533 0.7199 0.885 0.727 0.1057 0 0 0.8284 0.7116 0.4866 0 )
  let block-val-kat-2 ( list 0 0.5443 0.4624 0.5383 0.6492 0 0 0.5012 0.4789 0.4938 0.6435 0 0 0.5225 0.5076 0.4741 0.5791 0.5663 0 0.9995 0.6703 0.6761 0.8228 0.6815 0.871 0.9514 0.8679 0.5161 0.6802 0.7584 0.9344 0.9618 1 0.3191 0.2393 0.2157 0.815 0.813 0.9247 0.5042 0.212 0.1856 0.769 0.7131 0.711 0.2889 0.1534 0.118 0.7809 0.5221 0.5176 0.374 0.244 0.1908 0 0 0.4939 0.4498 0.5051 0 )
  let block-val-kat-3 ( list 0 0.5465 0.4406 0.4862 0.4404 0 0 0.4853 0.3784 0.4794 0.5802 0 0 0.567 0.4035 0.4414 0.5272 0.6176 0 1 0.5381 0.382 0.6403 0.6813 0.8465 0.9411 0.8096 0.5441 0.7883 0.7735 0.9179 0.8931 0.9637 0.3427 0.3232 0.3117 0.9205 0.8392 0.8785 0.4564 0.3102 0.305 0.8802 0.7424 0.6202 0.322 0.2576 0.2078 0.9401 0.5889 0.5387 0.4164 0.2382 0.1777 0 0 0.5698 0.4275 0.4118 0 )
  let block-val-kat-4 ( list 0 0.4315 0.3406 0.2776 0.5347 0 0 0.5089 0.3525 0.259 0.3855 0 0 0.6569 0.4462 0.3708 0.4415 0.7523 0 0.8659 0.5877 0.4377 0.6658 0.8918 0.8696 0.8685 0.9511 0.8375 0.7861 0.8468 0.9827 0.9134 0.932 0.6272 0.6809 0.5617 0.9968 1 0.9454 0.8803 0.6762 0.4845 0.9085 0.9043 0.8077 0.6977 0.5416 0.4581 0.8676 0.7344 0.7909 0.7304 0.64 0.6211 0 0 0.8921 0.8097 0.8367 0 )
  if ( area = "Homogeneous" ) [ set block-values block-val-homo ]
  if ( area = "Kattegat" and quarter = 1 ) [ set block-values block-val-kat-1 ]
  if ( area = "Kattegat" and quarter = 2 ) [ set block-values block-val-kat-2 ]
  if ( area = "Kattegat" and quarter = 3 ) [ set block-values block-val-kat-3 ]
  if ( area = "Kattegat" and quarter = 4 ) [ set block-values block-val-kat-4 ]
end


to landsc-display  ; Updates the display variable
  no-display
  if (disp-var = "bathymetry") [
    let min-bathymetry gis:minimum-of bathy-data
    let max-bathymetry gis:maximum-of bathy-data
    if (max-bathymetry - min-bathymetry = 0) [ set max-bathymetry max-bathymetry + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [ ; note the use of the "<= 0 or >= 0" technique to filter out
      ; "not a number" values, as discussed in the documentation.
      set pcolor 39
      if (bathymetry <= 0) or (bathymetry >= 0)
      [ set pcolor scale-color blue bathymetry max-bathymetry (min-bathymetry + 3) ]
    ]
  ]

  if (disp-var = "maxent-level") [ ; potential rate of food increase, from quarterly MAXENT preditions
    let min-food-gr gis:minimum-of maxent-level-data
    let max-food-gr gis:maximum-of maxent-level-data
    if (max-food-gr - min-food-gr = 0) [ set max-food-gr max-food-gr + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [
      set pcolor 39
      if maxent-level >= 0
      ;      [ set pcolor scale-color green maxent-level  ((log max-food-gr 10) + 1) ((log (min-food-gr + 0.1) 10) + 0.6)  ]
      [ set pcolor scale-color green maxent-level  1 -0.2 ] ; should work for probabilities
    ]
  ]

  if (disp-var = "food-prob") [ ; randomly distributed food patches
    let min-food-prob gis:minimum-of food-prob-data
    let max-food-prob gis:maximum-of food-prob-data
    ask patches
    [
      set pcolor 39
      if (food-prob <= 0) or (food-prob >= 0)
      [ set pcolor scale-color violet food-prob max-food-prob min-food-prob ]
    ]
  ]

  if (disp-var = "food-level") [ ; actual amount of food in patches
    let min-food-level 0
    let max-food-level maxU
    ask patches
    [
      set pcolor 39
      if maxent-level >= 0 [
        if ( food-level = 0 ) [ set pcolor white ]
        if ( food-level > 0 and food-level <= 0.02 * maxU ) [ set pcolor 45 ]
        if ( food-level > 0.02 and food-level <= 0.2 * maxU ) [ set pcolor 27 ]
        if ( food-level > 0.2 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 66 ]
        if ( food-level > 0.5 * maxU and food-level <= 0.9 * maxU ) [ set pcolor 63 ]
        if ( food-level > 0.9 * maxU and food-level < 1.1 * maxU ) [ set pcolor 61 ]
        if ( food-level > 1.1 * maxU ) [ set pcolor 61 ]
      ]
    ]
  ]
  if ( disp-var = "blocks" ) [
    let min-block 1
    let max-block 60
    ask patches
    [
      set pcolor 39
      if (block <= 0) or (block >= 0)
      [ set pcolor scale-color red block max-block min-block ]
    ]
  ]


;;; ********************* ADDED BY CARA ********************** ;;;

    ; visualize temperature values
    if (disp-var = "temperature") [
    let min-temp gis:minimum-of temp-data
    let max-temp gis:maximum-of temp-data
    if (max-temp - min-temp = 0) [ set max-temp max-temp + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [
      set pcolor 39
      if (temp-w <= 0) or (temp-w >= 0) ; filter out "NaN" values
      [set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 11 temp-w 20 0.1]
    ]
  ]
    ; visualize salinity values
    if (disp-var = "salinity") [
    let min-salinity gis:minimum-of salinity-data
    let max-salinity gis:maximum-of salinity-data
    if (max-salinity - min-salinity = 0) [ set max-salinity max-salinity + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [
      set pcolor 39
      if (salinity-w <= 0) or (salinity-w >= 0) ; filter out "NaN" values
      [set pcolor scale-color blue salinity-w max-salinity (min-salinity - 2)]
    ]
  ]

  ;;; ********************* End: ADDED BY CARA ********************** ;;;

  display
end  ; end landsc-display


to landsc-upd-food
  ; maxent-level is patch specific, between 0 and 1 (MAXENT-based); food-growth-rate (rU) is global variable
  ask patches with [ food-prob > 0 and food-level < (maxU * maxent-level) ] [
    let f-lev food-level + ( food-growth-rate * food-level * ( 1 - food-level / (maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters) ) )
    if (abs (f-lev - food-level) > 0.001) [
      repeat 47 [ set f-lev f-lev + ( food-growth-rate * f-lev * ( 1 - f-lev / (maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters) ) ) ]
    ]  ; If the food level is really low, let food grow 48 times -- like growing every half-hour step, only faster
    set food-level f-lev
    ; here maxent-level is MAXENT prediction and food-growth-rate is a universal calibrated variable
  ]

  ;if year >= 2 [ask patches with [ food-prob > 0 and food-level > (maxU * maxent-level * 1.1) ] [  set food-level food-level - (0.05) ]]

  ; update patch color - CHANGED BY CARA: TURNED OFF FOR PERFORMANCE
  ;  ask patches with [ food-prob > 0] [  if maxent-level >= 0 [
       ; if ( food-level = 0 ) [ set pcolor white ]
       ; if ( food-level > 0 and food-level <= 0.02 * maxU ) [ set pcolor 45 ]
       ; if ( food-level > 0.02 and food-level <= 0.2 * maxU ) [ set pcolor 27 ]
       ; if ( food-level > 0.2 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 66 ]
       ; if ( food-level > 0.5 * maxU and food-level <= 0.9 * maxU ) [ set pcolor 63 ]
       ; if ( food-level > 0.9 * maxU and food-level < 1.1 * maxU ) [ set pcolor 61 ]
       ; if ( food-level > 1.1 * maxU ) [ set pcolor 61 ]
  ;]]
end

to landsc-upd-maxent-level-map
  if ( area = "Kattegat" ) [
    let yyy word "quarter" quarter
    set yyy word yyy ".asc"
    set maxent-level-data gis:load-dataset word path yyy
    gis:apply-raster maxent-level-data maxent-level
    ask patches [ ifelse not ( maxent-level >= 0 ) [ set maxent-level -1 ] [ set maxent-level maxent-level ] ]  ; ### new check in PorpoiseModel version 1.1
    ask patches [ if food-prob > 0 and food-level < 0.01 [ set food-level 0.01 ] ]; Updated to ensure no zero food values
   ; landsc-display
  ]
end


; Porpoise variables and methods

to porps-setup-ref
  ; reference porpoises (with who = 0) -- deployment information
  if (area = "Kattegat") [
;    set deploy-x ( 619399.6 - xllcorner ) / 400  ; start-position
;    set deploy-y ( 6358048 - yllcorner ) / 400
    ; E of Sealand
    set deploy-x 527
    set deploy-y 397
    ; N Kattegat:
    set deploy-x 218
    set deploy-y 955
  ]
  if (area = "Homogeneous") [
    set deploy-x ( 619399.6 - xllcorner ) / 400  ; start-position, in pixels
    set deploy-y ( 6148048 - yllcorner ) / 400
  ]
  setxy deploy-x deploy-y
  set color red
  print "Porps deployed (red dot = porp 0)"
end

to porps-setup
  set memory-max 120
  set inertia-const 0.001
  set corr-logmov 0.94
  set corr-angle 0.26
  set vt list 0 0
  set ref-mem-strength-list-fixed ( list 0.999 0.9989 0.9988 0.9987 0.9985  0.9984  0.9982  0.9981  0.9979  0.9976  0.9974  0.9972  0.9969  0.9966  0.9962  0.9958  0.9954  0.995  0.9945  0.9939  0.9933  0.9926  0.9919  0.9911  0.9902  0.9893  0.9882  0.987  0.9858  0.9843  0.9828  0.9811  0.9793  0.9772  0.975  0.9726  0.9699  0.967  0.9638  0.9603  0.9565  0.9523  0.9478  0.9428  0.9375  0.9316  0.9252  0.9183  0.9108  0.9027  0.8939  0.8844  0.8742  0.8632  0.8514  0.8387  0.8252  0.8108  0.7954  0.7792  0.7619  0.7438  0.7248  0.7048  0.684  0.6624  0.64  0.617  0.5934  0.5692  0.5447  0.5199  0.4949  0.4699  0.445  0.4203  0.396  0.3721  0.3487  0.326  0.304  0.2828  0.2626  0.2432  0.2248  0.2074  0.1909  0.1755  0.161  0.1475  0.1349  0.1233  0.1125  0.1025  0.0933  0.0848  0.0771  0.0699  0.0634  0.0575  0.0521  0.0471  0.0426  0.0386  0.0349  0.0315  0.0284  0.0257  0.0232  0.0209  0.0189  0.017  0.0153  0.0138  0.0125  0.0112  0.0101  0.0091  0.0082  0.0074 )
  set work-mem-strength-list-fixed ( list 0.9990  0.9988  0.9986  0.9983  0.9979  0.9975  0.9970  0.9964  0.9957  0.9949  0.9938  0.9926  0.9911  0.9894  0.9873  0.9848  0.9818  0.9782  0.9739  0.9689  0.9628  0.9557  0.9472  0.9372  0.9254  0.9116  0.8955  0.8768  0.8552  0.8304  0.8022  0.7705  0.7351  0.6962  0.6539  0.6086  0.5610  0.5117  0.4617  0.4120  0.3636  0.3173  0.2740  0.2342  0.1983  0.1665  0.1388  0.1149  0.0945  0.0774  0.0631  0.0513  0.0416  0.0336  0.0271  0.0218  0.0176  0.0141  0.0113  0.0091  0.0073  0.0058  0.0047  0.0037  0.0030  0.0024  0.0019  0.0015  0.0012  0.0010  0.0008  0.0006  0.0005  0.0004  0.0003  0.0003  0.0002  0.0002  0.0001  0.0001  0.0001  0.0001  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 )
  let age-dist-list (list 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  5  5  5  5  5  5  5  5  5  5  5  5  5  5  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  7  8  8  8  8  8  8  8  8  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 10 10 10 11 11 11 11 11 12 12 12 12 12 12 12 13 13 13 13 14 14 14 14 15 15 15 15 18 18 19 19 21 22) ;stranded + bycaught animals in Lockyer & Kinze 2003 (NAMMCO)
  ;let birth-rate (list 0 0 0.136 0.417 0.818 0.714 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833)
  set max-movlgt 0
  set mean-movlgt 0
  set use-exp-food-val? false
  set CRW-contrib -9999
  ; set MR-contrib -9999
  set min-depth 1               ; water depth where porps very rarely swim (according to Jakob)
  set turn-right 1
  set movlgt-list list (0) (0)  ; two inputs required...  list used for histogramming
  set movlgt-list remove-item 0 movlgt-list
  set movlgt-list remove-item 0 movlgt-list
  set angle-list list (0) (0)  ; two inputs required...  list used for histogramming
  set list-of-dead-age [ ]
  set list-of-dead-day [ ]
    ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
  set list-of-dead-age-calves [ ]
  set list-of-dead-day-calves [ ]
   ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;
  set age-distrib [ ]
  set angle-list remove-item 0 angle-list
  set angle-list remove-item 0 angle-list
  create-porps n-porps
  ask porps [
      set age item (random (length age-dist-list)) age-dist-list
      set mating-day 0  ; set in yearly-tasks
      set ds-mating -99
      set dsg-birth -99
      set wean-scale-fact 1
      ifelse ( (age + 0.5) > age-of-maturity )[ set pregnancy-status 2 ][ set pregnancy-status 0 ] ; 0.5 added to correct for Jan 1st start date and integer age class
      if ( pregnancy-status = 2 and random-float 1 <= pregnancy-rate ) [  ;  become pregnant with prob. taken from Sørensen & Kinze 1994
          set pregnancy-status 1
          set ds-mating 360 - round( random-normal (7.5 * 360 / 12) 20 )

       ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
          set mass-f random-normal 0.74 0.19                               ; set initial fetal mass at approx 0.74 kg, corresponding to 135 ± 20 days into pregnancy using the relationship from Lockyer and Kinze 2003
       ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

      ]
      set with-lact-calf? false
      set disp-type 0
      set disp-target list 0 0
      set prev-logmov 0.8 ; unit = patchsize, i.e. times 100 m
      set prev-angle 10
      let jj 0
      let in-water? false
      while [ jj < 100 and in-water? = false ] [
        set deploy-x random-xcor
        set deploy-y random-ycor
        setxy deploy-x deploy-y
        if ( [ bathymetry ] of patch-here > 1 ) [ set in-water? true ]
        set jj jj + 1
      ]
      set enough-water-ahead? true
      let pos list deploy-x deploy-y
      set pos-list ( list pos )
      let tmp (list 0 0)
      set pos-list-daily ( list tmp tmp tmp tmp tmp tmp tmp tmp tmp tmp )
      ; set ref-mem-strength-list [ ]
      ; set work-mem-strength-list [ ]
      ; set work-mem-updated false
      set VE-total 0
      set deter-vt list 0 0
      set deter-strength 0
      set deter-time-left 0
      set stored-util-list  [ ]
      set color orange
      set shape "circle"
      set size 4
      if (who = 0) [ porps-setup-ref ]
      set pen-size 0.05

    ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
    porps-setup-params   ; porpoises then set their state variable values
    ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;
    ]

  carefully [
    set movlgt-list lput ( 10 ^ ( [ prev-logmov ] of porp 0 ) ) movlgt-list
  ]
  [
    set movlgt-list lput 1 movlgt-list
  ]
  if ( not (is-turtle? ( porp 0 )) ) [ porps-setup ]    ; strange that I have to do this... the porpoise isn't always created in the first go

;  landsc-display ; update displayed amount of food etc
end

to porps-upd-deter
  ask porps [
    if (deter-time-left = 0) [
      set deter-strength 0
      set deter-vt replace-item 0 deter-vt 0
      set deter-vt replace-item 1 deter-vt 0
    ]
    if (deter-time-left > 0) [
      ; reduce deterrence by half every time step
      set deter-time-left deter-time-left - 1
      set deter-strength deter-strength / 2  ;
      set deter-vt replace-item 0 deter-vt ( item 0 deter-vt / 2 )
      set deter-vt replace-item 1 deter-vt ( item 1 deter-vt / 2 )
    ]
  ]
end ; end porps-upd-deter

to porp-check-depth
  ; Check that there is enough water at all steplengths ahead, set enough-water-ahead? to false if < min-depth
  set enough-water-ahead? true
  let pres-mov ( 10 ^ pres-logmov )                                                                           ; because pres-logmov may have changed in the porp-avoid-land procedure
  let dd ceiling ( pres-mov / 0.1 )                                                                           ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [ bathymetry ] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ?1 -> ?1 > 0 ] depth-list )) ) [                               ; i.e. if some depths on the list aren't > 0
    set enough-water-ahead? false
  ]
end


to porp-avoid-land
  ; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
  ; Don't do the turning here, but change angle to be turned in porp-std-move or porp-markov-mov.
  ; Note that the emergency procedure "avoid-beh 5" is found in porp-std-move
  let rand-ang random 10
  let avoid-beh 0
  let pres-mov ( 10 ^ pres-logmov )
  let bath-l [ bathymetry ] of patch-left-and-ahead (40 + rand-ang) pres-mov
  let bath-r [ bathymetry ] of patch-right-and-ahead (40 + rand-ang) pres-mov
  ; alternative kinds of evasive behaviour:
  if ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 1  ; evasive behaviour type 1
    ifelse ( bath-r >= min-depth and bath-l >= min-depth )
      [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
      [ ifelse ( bath-r >= min-depth )
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
  ]
  ; else try turning more aprubtly ( = 70 deg )
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 2  ; evasive behaviour type 2
    set bath-l [ bathymetry ] of patch-left-and-ahead (70 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (70 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth )
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth )
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
    ]
  ]
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 3  ; evasive behaviour type 3
    set bath-l [ bathymetry ] of patch-left-and-ahead (120 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (120 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth )
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth )
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
    ]
  ]
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    ; if everything else fails, turn around
    set avoid-beh 4  ; evasive behaviour type 4
    let j 0
    if deep-water? = false [ porp-check-depth ]                     ; CHANGED BY CARA
    while [ not enough-water-ahead? and j < length pos-list ] [
      facexy (item 0 (item j pos-list)) (item 1 (item j pos-list))  ; each item on pos-list contains a list with a x and a y-coordinate
      setxy (item 0 (item j pos-list)) (item 1 (item j pos-list))
      set j j + 1
      if deep-water? = false [ porp-check-depth ]    ; CHANGED BY CARA
      if (j = 20) [ set enough-water-ahead? true ]
    ]
  ]
  if ( debug = 1 ) [
    let tmp-list list ("beh =") avoid-beh
    set tmp-list lput ("; tck =") tmp-list
    set tmp-list lput time-step tmp-list
    write tmp-list
    let tmp7 word "; " round pres-angle
    print word tmp7 " degr."
  ]
end ; end porp-avoid-land


to porp-markov-move
  ; Movements based on dead-reckoning data -- first calc distance, then turning angle
  set pres-logmov ( 0.5 + random-normal 0 0.25 )
  let pres-mov ( 10 ^ pres-logmov )
  set pres-angle random-normal 0 40
  if ( abs pres-angle > 60 ) [ set pres-angle (1 + random-float 0.5) * pres-angle ]  ; make angle dist more leptokurtic
  right pres-angle
  ;
  ; Turn to avoid swimming on land if necessary:
  ; ( section copied in from porp-std-move)
  let goto-avoid-land? false
  let dd 0

  ifelse deep-water? = false [    ; CHANGED BY CARA
  set dd ceiling ( pres-mov / 0.25 )  ; number of 25-m steps to check water depth at
  set goto-avoid-land? false
  if (not ( [ bathymetry ] of patch-ahead pres-mov >= min-depth ) ) [ set goto-avoid-land? true ]
  repeat dd [
    if ( not ( [ bathymetry ] of patch-ahead ( dd * 0.25 ) >= min-depth ) ) [ set goto-avoid-land? true ]  ; must use "not >= " rather than " < " for catching NaN's
    set dd dd - 1
    ]
  ]
  [ set  goto-avoid-land? false ]

  if ( goto-avoid-land? ) [ porp-avoid-land ]

  if deep-water? = false [   ; CHANGED BY CARA
  set pres-mov ( 10 ^ pres-logmov )  ; because pres-logmov may have changed in the porp-avoid-land procedure
;  let ready-to-move true
  ; test again:
  set dd ceiling ( pres-mov / 0.1 )  ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [bathymetry] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ?1 -> ?1 > 0 ] depth-list )) ) [ ; i.e. if some items on the list aren't < 0
    uphill bathymetry
    if ( debug = 1 ) [
      show word "Tick = " time-step
      show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here)
    ]
    ]
  ]
  ;
  ; move
;  if (ready-to-move) [
  fd pres-mov / 4  ; divide by four when cell size is 400 x 400 m
;     ]
  ; Remember current moves for the next iteration
  set pres-logmov log pres-mov 10
  set prev-angle pres-angle
  set prev-logmov pres-logmov
end


to porp-ref-mem-turn
  ; Move towards places visited previously if food was found there and they aren't too far away or forgotten.
  let bb ( [food-level] of patch-here )  ; Stationary food species. The stored intrisic patch utility for t=0. Initially it is either 0, 1,
  ; or -9999, but grows logistically after food is eaten
  set vt list -9999 -9999  ; init to catch errors
  if not ( abs(bb) >= 0 ) [
    ; There are errors in food availability -- sometimes Na is calculated even though depth is > 0. Catch error here
    set bb 0
    if ( debug = 4 ) [
      print "Replaced NaN food value with 0"
      print patch-here
    ]
  ]
  set stored-util-list fput bb stored-util-list

  ; Update reference memory strength for past locations
  ; Replaced by "ref-mem-strength-list-fixed"
  let max-mem 0.999
  ; set ref-mem-strength-list fput max-mem ref-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  ; let ii 1
  ; while [ ii < length ref-mem-strength-list ] [
  ;   let MRPt item (ii - 1) ref-mem-strength-list  ; (reference memory for patch p at time t)
  ;   let reduced-mem MRPt - ref-mem-decay * (1 - MRPt) * MRPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
  ;   set ref-mem-strength-list replace-item ii ref-mem-strength-list reduced-mem
  ;   set ii ii + 1
  ; ]

  ; Set patch value for each past location -- perceived patch utility (= reference memory x intrinsic patch utility (stuff eaten)), divided by DIRECT distance
  let perceived-util-list [ ]
  set perceived-util-list lput (item 0 stored-util-list * max-mem) perceived-util-list
  let tmp list (0) (0)
  let attr-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (this place) has length 0, the others have length 1
  set attr-vector-list lput tmp attr-vector-list
  let one-attr-vector [ ]
  let vector-lgt 0
  let ii 1
  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    if (item ii stored-util-list = 0) [  ; save time by skipping dist measure when there is no food -- changed 101214
      set perceived-util-list lput 0 perceived-util-list
      set one-attr-vector list 0 0
      set attr-vector-list lput one-attr-vector attr-vector-list
    ]
    if (not ( item ii stored-util-list = 0 )) [  ; save time by skipping dist measure when there is no food -- changed 101214
      set dist-to-foodpos (distancexy ( item 0 (item ii pos-list) ) ( item 1 (item ii pos-list) ))
      ifelse (dist-to-foodpos < 1E-20 )
        [ set perceived-util-list lput 9999 perceived-util-list ]      ; arbitrary large value for close dist
        [ set perceived-util-list lput ( (item ii stored-util-list) * (item ii ref-mem-strength-list-fixed) / dist-to-foodpos ) perceived-util-list ]
      ; = utility * memory / distance
      ; Create attraction vectors; unit-vectors pointing towards the patches in memory
      set one-attr-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
      ; make sure that it works with wrapping landscapes:
      if ( item 0 one-attr-vector > world-width / 2 ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector - world-width ) ]
      if ( item 0 one-attr-vector < (- world-width / 2 ) ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector + world-width ) ]
      if ( item 1 one-attr-vector > world-height / 2 ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector - world-height ) ]
      if ( item 1 one-attr-vector < (- world-height / 2 ) ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector + world-height ) ]
      set vector-lgt  sqrt( item 0 one-attr-vector * item 0 one-attr-vector + item 1 one-attr-vector * item 1 one-attr-vector )
      if vector-lgt = 0 [
        if ( debug = 4 ) [
          show word "attr-vector-lgt = " vector-lgt
          print "skipping to next porp"
        ]
        stop
      ]
      set one-attr-vector replace-item 0 one-attr-vector ((item 0 one-attr-vector) / vector-lgt)
      set one-attr-vector replace-item 1 one-attr-vector ((item 1 one-attr-vector) / vector-lgt)
      set attr-vector-list lput one-attr-vector attr-vector-list
    ]
    set ii ii + 1
  ]
  ; Calculate resultant attraction vector vt as sum of products of individual values and attraction vectors (eqn 5). May have length != 1
  set ii 1  ; no attraction to current pos (t=0)
  let vt-x 0
  let vt-y 0
  while [ ii < length pos-list ] [
    set vt-x vt-x + item ii perceived-util-list * item 0 ( item ii attr-vector-list )
    set vt-y vt-y + item ii perceived-util-list * item 1 ( item ii attr-vector-list )
    set ii ii + 1
  ]
  if ( debug = 4 ) [
    type word "Food here: " bb
    type ",  Attr.v: "
    let attr-vect list vt-x (",")
    set attr-vect lput vt-y attr-vect
    print attr-vect
    if (not ( abs(vt-x) >= 0)) [  ; catch missing values
      write "Perc.util: "
      print perceived-util-list
    ]
  ]
  set vt list vt-x vt-y

  ; Remove items in distant past to increase execution speed
  ; if ( length ref-mem-strength-list-fixed > memory-max ) [ set ref-mem-strength-list remove-item memory-max ref-mem-strength-list ]
  if ( length stored-util-list > memory-max ) [ set stored-util-list remove-item memory-max stored-util-list ]
end  ; end porp-ref-mem-turn


to porp-work-mem-turn
   ; Influences direction moved in std-move through vector 'vt'
   ; This procedure MUST be called after porp-ref-mem-turn, as it uses the stored-util-list (experienced food at patch) calculated there,
   ; and because it updates the vt vector which is later used in std.move (adds to the previously calculated vt).

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  ; set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009


;  print work-mem-strength-list ; ### TMP
  let ii 1
;  if ( work-mem-updated = false ) [
;    while [ ii < length work-mem-strength-list ] [
;      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
;      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
;      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
;      set ii ii + 1
;    ]
;  ]
;  set work-mem-updated true
  ; Presently no need to multiply with stored-util-list to get perceived utility -- the animal knows that all food is eaten in the patch.
;  print work-mem-strength-list ; ### TMP


  let tmp list (0) (0)
  let deter-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (0; this place) has length 0, the others have length 1
  set deter-vector-list lput tmp deter-vector-list
;  show word "deter-vector-list: " deter-vector-list
  let one-deter-vector [ ]
  let vector-lgt 0
  set ii 1
;  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    ; Create deterrence vectors; unit-vectors pointing towards the patches in memory
    set one-deter-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-deter-vector > world-width / 2 ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector - world-width ) ]
    if ( item 0 one-deter-vector < (- world-width / 2 ) ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector + world-width ) ]
    if ( item 1 one-deter-vector > world-height / 2 ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector - world-height ) ]
    if ( item 1 one-deter-vector < (- world-height / 2 ) ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-deter-vector * item 0 one-deter-vector + item 1 one-deter-vector * item 1 one-deter-vector )
    if vector-lgt = 0 [
      if ( debug = 5 ) [
        show word "deter-vector-lgt = " vector-lgt
        print "Haven't moved, skipping to next porp"
      ]
      stop
    ]
    set one-deter-vector replace-item 0 one-deter-vector ((item 0 one-deter-vector) / vector-lgt)
    set one-deter-vector replace-item 1 one-deter-vector ((item 1 one-deter-vector) / vector-lgt)
    set deter-vector-list lput one-deter-vector deter-vector-list
    set ii ii + 1
  ]

  ; Calculate resultant deterrence vector vtd as sum of products of individual values and deterrence vectors
  set ii 1  ; no deterrence from current pos (t=0)
  let vtd-x 0
  let vtd-y 0
  while [ ii < length pos-list ] [
    set vtd-x vtd-x + inertia-const * item ii work-mem-strength-list-fixed * item 0 ( item ii deter-vector-list )
    set vtd-y vtd-y + inertia-const * item ii work-mem-strength-list-fixed * item 1 ( item ii deter-vector-list )
    set ii ii + 1
  ]

  if ( debug = 5 ) [
;    print word "work-mem: " work-mem-strength-list
    type "Deter.v: "
    let deter-vect list vtd-x (",")
    set deter-vect lput vtd-y deter-vect
    print deter-vect
    if (length pos-list > 1) [ print word "pos. before: " item 1 pos-list ]
    print word "pos. now: " item 0 pos-list
    print ""
    ; Checked -- works, at least with length pos-list = 2
  ]

  ; vtd points towards the previous position, must be subtracted from vt
  set vt replace-item 0 vt ( item 0 vt - vtd-x )
  set vt replace-item 1 vt ( item 1 vt - vtd-y )

; Remove items in distant past to increase execution speed
;  if ( length work-mem-strength-list > memory-max ) [ set work-mem-strength-list remove-item memory-max work-mem-strength-list ]
end ; end porp-work-mem-turn


to porp-get-exp-food-val
  ; Calculate the expaected value (VE-total) of the food to be found in the future based on food found in recent positions x the working memory
  ; Uses the values of the patches in "stored-util-list", calculated in porp-ref-mem-turn

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
;  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
;  let ii 1
;  if ( work-mem-updated = false ) [  ; list may have been updated in porp-work-mem-turn
;    while [ ii < length work-mem-strength-list ] [
;      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
;      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
;      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
;      set ii ii + 1
;    ]
;  ]
;  set work-mem-updated true

  let ii 1
  set VE-total 0
  let max-i min list ( length work-mem-strength-list-fixed ) ( length stored-util-list )
  while [ ii < max-i ] [
    set VE-total VE-total + item (ii - 1) work-mem-strength-list-fixed * item (ii - 1) stored-util-list
    set ii ii + 1
  ]

  if ( debug = 5 ) [
    print ""
    ; print word "stored-util-list: " stored-util-list
    ; print word "work-mem-strength-list: " work-mem-strength-list
    let ttt (list "VE-total for porp " who ":")
    show word ttt VE-total
  ]
end ; end porp-get-exp-food-val


to porp-std-move
; Movements based on dead-reckoning data:
  ; global vars: corr-logmov 0.94 and corr-angle 0.26

  ; ### turning angle
  let prev-mov 10 ^ prev-logmov
  let pres-heading heading
  set pres-angle 999
  let j 1
  let tmp-angle 0
  ifelse ( prev-angle < 0 ) [ set tmp-angle prev-angle - 24 ] [ set tmp-angle prev-angle + 24 ]  ; for increasing mean turning angle
  while [ abs (pres-angle) > 180 ]  [
    let R2 random-normal 0 38
    set pres-angle ( tmp-angle * (- corr-angle) + R2 )                  ; Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water.
    set j j + 1
    if (j = 200) [
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [
        print word "exiting loop 1, ang=" pres-angle
      ]
    ]
  ]
  let sign 0
  ifelse pres-angle < 0 [ set sign -1 ] [ set sign 1 ]
  set pres-angle abs pres-angle ; add the sign again later
  ; Make angle decrease linearly with mov-dist
  let go-on? true
  set j 1
  let rnd 0
  while [ go-on? ]  [
    set rnd random-normal 96 28      ; draws the number to be added to pres-angle
    if ( prev-mov <= 5.50 ) [
      set pres-angle pres-angle + rnd - ( rnd * prev-mov / 5.50 )
    ]
    if ( pres-angle < 180 ) [ set go-on? false ]  ; remember that turning angle is unsigned here
    set j j + 1
    if (j = 200) [
      set pres-angle ( random 20 + 90 )
      set go-on? false
      if ( debug = 3 ) [
        print word "exiting loop 2, ang=" pres-angle
      ]
    ]
  ]
  ; if ( abs pres-angle > 55 and abs pres-angle < 180 ) [ set pres-angle (1 - random-float 0.32) * pres-angle ]  ; make turning angle dist more leptokurtic
  set pres-angle pres-angle * sign
  let angle-before-avoid-land pres-angle ; for printing later using debug 2
  right pres-angle
  let angle-turned-right pres-angle ; for updating prev-angle at end of porp-std-move
  set pres-angle 0

  ; ### distance
  set pres-logmov 999
  let porp-max-dist 1.18                                                                                        ; log10 ( max distance a porpoise can travel per half-hour )
  set j 1
  while [ pres-logmov > porp-max-dist ] [
    let R1 random-normal 0.42 0.48
    set pres-logmov ( corr-logmov * prev-logmov + R1 )
    set j j + 1
    if (j = 200) [
      if (pres-angle = 0) [set pres-angle pres-angle + 0.00001]
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [
        print word "exiting loop 3, ang=" pres-angle
      ]

    ]
  ]
  let pres-mov ( 10 ^ pres-logmov )                                                                             ; This is what is plotted in the histogram
  ;
  ; Turn to avoid swimming on land if necessary:
  ifelse deep-water? = false [ set enough-water-ahead? false ] [ set enough-water-ahead? true ]    ; CHANGED BY CARA
  let count-i 0
  while [ not enough-water-ahead? ] [
    porp-check-depth
    if (not enough-water-ahead?) [ porp-avoid-land ]
    set pres-mov ( 10 ^ pres-logmov )                                                      ; because pres-logmov may have changed in the porp-avoid-land procedure
    right pres-angle                                                                       ; angle to turn -- pres-angle -- is changed in porp-avoid-land
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-angle 0
    set count-i count-i + 1
    if (count-i = 100) [
      set enough-water-ahead? true
      if ( debug = 1 ) [
        print "caught water-ahead loop"
      ]
    ]
  ]
  ; test depth again, avoid-beh = 5:
  ifelse deep-water? = false [     ; CHANGED BY CARA
    set enough-water-ahead? false
    porp-check-depth
  ]
  [ set enough-water-ahead? true ]

  if (not enough-water-ahead?) [
    let prev-heading heading
    let p max-one-of neighbors [ bathymetry ]
    carefully [ face p ] [ ]
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-mov 1                                                                                            ; move 100 m towards deeper patch
    if ( debug = 1 ) [
      let tmp5 list "beh =  5 ; tck " time-step
      write tmp5
    ]
  ]
  ;
  ; slow down if turning sharply:
  if ( pres-mov > 10 and (abs angle-turned-right) > 90 ) [ set pres-mov  pres-mov / 5  ]
  if ( pres-mov > 7 and (abs angle-turned-right) > 50 ) [ set pres-mov  pres-mov / 2  ]
  ;
  ; Change direction if attracted / deterred by certain areas (model >= 2)
  let total-dx 0
  let total-dy 0
  if ( not use-exp-food-val? ) [
    set total-dx (dx * pres-mov) + (item 0 vt)                ; vt isn't used in porp-std-move till here
    set total-dy (dy * pres-mov) + (item 1 vt)                ; note that dx is change in x if taking ONE step forward
    facexy (xcor + total-dx) (ycor + total-dy)
  ]
  if ( use-exp-food-val? and model < 3  ) [
    set CRW-contrib inertia-const + pres-mov * VE-total       ; length of vector pointing in direction predicted by CRW (VE-total and pres-mov are porp variables)
    ; set MR-contrib sqrt ( (item 0 vt) * (item 0 vt) + (item 1 vt) * (item 1 vt) )     ; length of vector pointing in direction of remembered food
    set total-dx (dx * CRW-contrib) + (item 0 vt)
    set total-dy (dy * CRW-contrib) + (item 1 vt)
    facexy (xcor + total-dx) (ycor + total-dy)                ; really not needed, it already points that way
  ]
  ; deterrence behaviour -- get scared away from ships and wind turbines
  if ( use-exp-food-val? and model >= 3 ) [
    set CRW-contrib inertia-const + pres-mov * VE-total
    set total-dx (dx * CRW-contrib) + (item 0 vt) + ((item 0 deter-vt) * deterrence-coeff)
    set total-dy (dy * CRW-contrib) + (item 1 vt) + ((item 1 deter-vt) * deterrence-coeff)
    facexy (xcor + total-dx) (ycor + total-dy)
  ]
  ; for testing - CARA::::
  ;let tot-dx-CRW (dx * CRW-contrib) + (item 0 vt)
  ;let tot-dy-CRW (dy * CRW-contrib) + (item 1 vt)
  ;set tot-vt-CRW sqrt ((tot-dx-CRW * tot-dx-CRW) + (tot-dy-CRW * tot-dy-CRW))

  ;let tot-dx-deter ((item 0 deter-vt) * deterrence-coeff)
  ;let tot-dy-deter ((item 1 deter-vt) * deterrence-coeff)
  ;set tot-vt-deter sqrt ((tot-dx-deter * tot-dx-deter) + (tot-dy-deter * tot-dy-deter))

  ; Store turn for calc of turning angle in next step:
  ; let total-turn heading - pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move
  let total-turn subtract-headings heading pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move

  ;
  ; Move:
  ; In the population model all movement lengths are still calculated in 100 m steps, but the cell size has increased from 100 m to 400 m
  ; The step should therefore be divided by 4
  fd pres-mov / 4  ; movement length isn't affected by presence of food
  ;
  if ( debug = 2 ) [
    if ( time-step = 0 ) [ print "dist angle-before-avoid-land angle-turned-right x y" ]
    let tmp-var2 (round ( (10 ^ prev-logmov) * 100) / 100)    ; THIS IS IMPORTANT -- the porp turns before it moves, so turning angle is affected by previous moving dist
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-before-avoid-land
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-turned-right
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( xcor * 100 + xllcorner )
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( ycor * 100 + yllcorner )
    print tmp-var2
  ]
  if ( debug = 5 ) [
    print word "CRW-contrib: " list ( (dx * (inertia-const + pres-mov * VE-total)) ) ( (dy * (inertia-const + pres-mov * VE-total)) )
    print word "MR-contrib: " vt
    print word "dx, dy (after): " list (total-dx) (total-dy)
    print word "heading (after): " heading
    print word "total-turn: " heading
  ]
  ;
  ; Remember current moves for the next iteration
  ; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
  ; set prev-angle angle-turned-right  ; so the additional turn due to attraction to food does not influence turning angle in next step
  set prev-angle total-turn  ; so the additional turn due to attraction to food DOES influence turning angle in next step
  set prev-logmov log pres-mov 10  ; total steplength, resulting from vt + pres-mov
  ;
  ; test depth one last time, avoid-beh = 6 - move back on same track:
  if (not ([ bathymetry ] of patch-here > 0) ) [
    let prev-heading heading
    if (length pos-list > 1 ) [ facexy (item 0 item 1 pos-list) (item 1 item 1 pos-list) ] ; debugging 110323
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    if (length pos-list > 1 ) [ setxy (item 0 item 1 pos-list) (item 1 item 1 pos-list) ]                             ; move 100 m towards deeper patch
    if ( debug = 1 ) [
      ; print "; -> "
      let tmp6 list "beh =  6 ; tck " time-step
      set tmp6 word tmp6 " ; "
      set tmp6 word tmp6 angle-turned-right
      print word tmp6 " degr."
    ]
  ]
  ; update position list:
  let pres-pos list xcor ycor
  set pos-list fput pres-pos pos-list
  if ( length pos-list > memory-max ) [ set pos-list remove-item memory-max pos-list ]
end ;  end porp-std-move


; REMOVED PORP-UPD-ENERGETIC-STATUS - CHANGED BY CARA

to porp-move
  if ( model = 0 ) [ porp-markov-move ]
  if ( model = 1 ) [ porp-std-move ]
  if ( model = 2 ) [
    set use-exp-food-val? true
    porp-ref-mem-turn       ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt' (global var)
    porp-get-exp-food-val   ; determines the tendency to move following CRW behaviour based on foraging success in recent past
    porp-std-move           ; this is where the porp moves forward
    ;porp-upd-energetic-status           ; food level increases in 'go' -- affect the landscape and energetic status of the porpoise - CHANGED BY CARA
  ]
  if ( model >= 3 ) [
    set use-exp-food-val? true
    porp-ref-mem-turn       ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt' (global var)
    porp-get-exp-food-val   ; determines the tendency to move following CRW behaviour based on foraging success in recent past
    porp-std-move           ; this is where the porp moves forward and responds to noise by turning away
    ;porp-upd-energetic-status  ; transform food to energy and spend energy based on step length. Food level in patches increases in 'go' - CHANGED BY CARA
                            ; mortality and pregnancy status is set in daily-tasks for models > 4
  ]
  if not ( [ bathymetry ] of patch-here > 0 ) [
    follow-me
    beep
    user-message "Error, no water"
  ]
end ; end porp-move


to porp-disp-target-select
  ; deciding where to disperse to based on knowledge of other blocks (each block is 100 x 100 cells = 40 000 x 40 000 m)
  let nbr 0
  let block-quality ( list 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
;  if (debug = 7 and (who = 0 or who = 1)) [ print block-quality ]
  while [nbr < length (block-quality) ] [
    let one-dx ( abs (xcor - item nbr block-centres-x) )   ; x dist from porp to block centres -- distancexy doesn't work here
    let one-dy ( abs (ycor - item nbr block-centres-y) )
    let block-dist round ( sqrt( one-dx * one-dx + one-dy * one-dy ) )
;    if (debug = 7 and (who = 0 or who = 1)) [
;      write (nbr + 1)
;      print word "block-dist: " block-dist
;    ]
    set block-quality replace-item nbr block-quality ( item nbr block-values )
    ; set quality to 0 for blocks that are too close or too far (divide by 0.4 to convert from km to cells)
    if (block-dist < min-dist-to-target / 0.4 ) [ set block-quality replace-item nbr block-quality 0 ]
    set nbr nbr + 1
  ]
;  if (debug = 7 and (who = 0 or who = 1)) [ print block-quality ]

  ; numbers of the blocks with highest quality (blocks numbered 0-59, best block first)
  let hi-q-blocks [ ]
  if (n-disp-targets = 4) [ set hi-q-blocks (list 0 0 0 0 ) ]
  if (n-disp-targets = 6) [ set hi-q-blocks (list 0 0 0 0 0 0 ) ]
  if (n-disp-targets = 8) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0) ]
  if (n-disp-targets = 10) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0) ]
  if (n-disp-targets = 12) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0 0 0) ]
  let block-quality-srt sort block-quality   ; quality of the blocks. Last element corresponds to first element in hi-q-blocks
  let i 0
  while [i < length(hi-q-blocks)] [
    let hi-qual-i item ((length block-quality-srt) - i - 1) block-quality-srt  ; block w highest value
    let j 0
    while [ j < length(block-values) ] [
      if (item j block-quality = hi-qual-i) [
        set hi-q-blocks replace-item i hi-q-blocks j  ; number of best blocks in left end
      ]
      set j j + 1
    ]
    set i i + 1
  ]

  ; select block at random from the twelve blocks with highest quality (where qual = mean.food / dist)
  ; Make sure that porps far north do not try to disperse west
  let the-nbr -9
  set the-nbr random (length hi-q-blocks)
  let sel-block item the-nbr hi-q-blocks
  if ( [ block ] of patch-here < 20 and ( sel-block = 30 or sel-block = 31 or sel-block = 36 or sel-block = 37 or sel-block = 43) ) [ set sel-block sel-block + 2 ]  ; find block slightly more to the east if currently far north
  set disp-target list ( item sel-block block-centres-x ) ( item sel-block block-centres-y )  ; coordinates of cell to move towards

  if ( debug = 7 and (who = 0 or who = 1) ) [
    print " "
    let tmp word "Disp from: " [ block ] of patch-here
    set tmp word tmp " -> "
    show word tmp sel-block
  ]
 end ; end porp-disp-target-select


to porp-disp-1  ; disperse to new region
  ; When the porp get closer to land than min-dist-to-land it tries to turn away, else it shifts to disp-2
  let min-dist-to-land 2000

  ; If porp is SE of Sealand, don't do directed dispersal
  if (xcor > 438 and ycor < 478 ) [
  ;if ([block] of patch-here = 35 or [block] of patch-here = 41 ) [
    set disp-type 2
    stop
  ]

  ; If porp is N of Djursland or in Little Belt, don't do directed dispersal
  if ([block] of patch-here = 14 or [block] of patch-here = 31 ) [
    set disp-type 2
    stop
  ]

  ; Go north if north of Fyn / Funen
  if ([block] of patch-here = 32 ) [
    set heading 0
    set disp-type 2
    fd 1
    stop
  ]

  if ( debug = 7 and who = 0 ) [
      ask ( porp 0 ) [ pen-down ]
  ]
  if ( debug = 7 and who = 1 ) [
      ask ( porp 1 ) [ pen-down ]
  ]
  ; Find distance to target block centre
  let the-dx (item 0 disp-target) - xcor
  let the-dy (item 1 disp-target) - ycor
  let the-dist round ( sqrt( the-dx * the-dx + the-dy * the-dy ) )
  facexy ( xcor + the-dx / 2 ) ( ycor + the-dy / 2 )
;  if ( debug = 7 and (who = 0 or who = 1) ) [ print word "Heading A: " heading ]

  ; adjust angle to swim towards deep areas
  let bathymetry-ahead ( list
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist)
  )

  let bathymetry-far-ahead ( list
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist * 8)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist * 8)
  )

  ; Turn up to 20 degr towards deepest water, provided that there is no land further away in that direction
  let i 0
  let good-heading? ( list 0 0 0 0 0 0 0 0 0 )
  while [ i < length good-heading? ] [
    ifelse ( item i bathymetry-far-ahead > 0 ) [ set good-heading? replace-item i good-heading? true ] [ set good-heading? replace-item i good-heading? false ]
    set i i + 1
  ]

  let angles (list -40 -30 -20 -10 0 10 20 30 40 )
  let bathymetry-choice -999
  set i 0
  let sel-angle 0

  set bathymetry-choice max bathymetry-ahead
  while [ i < length bathymetry-ahead ] [
    if ( item i bathymetry-ahead = bathymetry-choice ) [
      set sel-angle item i angles]
    set i i + 1
  ]

  set heading heading + sel-angle

;  if (debug = 7 and (who = 0 or who = 1) ) [
;    type word "Heading B: " heading
;    print word ", angle: " sel-angle
;  ]

  ; Turn to areas far from land if there is land ahead
  let disttocoast-ahead ( list
    -999 ; [ disttocoast ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist * 2)
    -999 ; [ disttocoast ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist * 2)
  )

  ; make sure that there is also water far away
  set i 0
  while [ i < length(disttocoast-ahead) ] [
    if ( not item i good-heading? ) [ set disttocoast-ahead replace-item i disttocoast-ahead -999 ]
    set i i + 1
  ]

  let disttocoast-choice -999
  let low-water-ahead not (item 1 bathymetry-ahead > min-disp-depth and item 2 bathymetry-ahead > min-disp-depth and item 3 bathymetry-ahead > min-disp-depth and item 4 bathymetry-ahead > min-disp-depth and item 5 bathymetry-ahead > min-disp-depth) ; that is, if some pos ahead are on land
  if ( low-water-ahead or [ disttocoast ] of patch-here < min-dist-to-land ) [
    set disttocoast-choice max disttocoast-ahead
    set i 0
    while [ i < length disttocoast-ahead ] [
      if ( item i disttocoast-ahead = disttocoast-choice ) [ set sel-angle item i angles ]
      set i i + 1
    ]
  ]
  set heading heading + sel-angle
;  if (debug = 7 and (who = 0 or who = 1) ) [
;    type word "Heading C: " heading
;    print word ", angle: " sel-angle
;  ]



  ; ### when to stop dispersing:  ###

  ; shift to dispersal away from here, along coast, if porp cannot move for one day
  if ( distancexy (item 0 (item 1 pos-list-daily)) (item 1 (item 1 pos-list-daily)) < 2 ) [
    set disp-type 2
    if (debug = 7 and who = 0 ) [
      print "Not mov 1 d, chg to disp-type 2 (porp 0)"
    ]
    if (debug = 7 and who = 1 ) [
      print "Not mov 1 d, chg to disp-type 2 (porp 1)"
    ]
  ]

  ; shift to dispersal away from here, along coast, if porp moves too little
  if ( distancexy (item 0 (item 8 pos-list-daily)) (item 1 (item 8 pos-list-daily)) < 6 ) [
    set disp-type 2
    if (debug = 7 and (who = 0 or  who = 1) ) [
      show "Mov <15 km in 8d, chg to disp-type 2 (porp 1)"
    ]
  ]

  ; Close to coast, chg. disp. mode
  if ( [ disttocoast ] of patch-here < min-dist-to-land ) [
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1 )) [
      show "TOO CLOSE to land, chg to disp-type 2"
    ]
  ]

  ; Close to target, chg. disp. mode
  if ( the-dist < 50 ) [   ; each block is 100 x 100 cells
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1) ) [
      show "Close to target, chg to disp-type 2"
    ]
  ]

  if ( not ( [ bathymetry ] of patch-ahead (mean-disp-dist * 4) > 0 and [ bathymetry ] of patch-ahead (mean-disp-dist * 3) > 0 and [ bathymetry ] of patch-ahead (mean-disp-dist * 2) > 0 ) ) [
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1) ) [
      show "LAND ahead, chg to disp-type 2 (porp 1)"
    ]
  ]

  if ( not enough-water-ahead? ) [
    set disp-type 0
    if (debug = 7 and who = 0 ) [
      print "NO WATER, stop dispersing (1) (porp 0)"
      pen-up
    ]
    if (debug = 7 and who = 1 ) [
      print "NO WATER, stop dispersing (1) (porp 1)"
      pen-up
    ]
  ]

  if ( disp-type = 1 ) [
    fd ( mean-disp-dist / 0.4 )
    if not ([bathymetry] of patch-here >= 0) [ ; catching rare errors
      fd (- mean-disp-dist / 0.4)
      set disp-type 0
    ]
  ;  set energy-level energy-level - 0.001 * e-use-per-km * mean-disp-dist / 0.4   ; distances are measured in number of cells here ; REMOVED BY CARA
  ]

end ; end porp-disp-1


to porp-disp-2  ; disperse along coast, away from prev position
  ; Disperse at least this dist from coast (same var name as in porp-disp-1)
  let min-dist-to-land 500

  if (debug = 7 and (who = 0) ) [
    set color blue
    ask ( porp 0 ) [ pen-down ]
  ]
  if (debug = 7 and (who = 1) ) [
    set color blue
    ask ( porp 1 ) [ pen-down ]
  ]

  ; Turn away from place visited 1 day ago
  facexy (item 0 (item 1 pos-list-daily)) (item 1 (item 1 pos-list-daily))
  rt 180

  ; adjust angle to swim at const dist from land
  let disttocoast-ahead ( list
    [ disttocoast ] of patch-at-heading-and-distance (heading - 80) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 70) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 60) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 50) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 40) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 30) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 20) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 10) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 10) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 20) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 30) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 40) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 50) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 60) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 70) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 80) mean-disp-dist
  )


;  let low-dtc-ahead not (item 1 disttocoast-ahead > 500 and item 2 disttocoast-ahead > 500 and item 3 disttocoast-ahead > 500 and item 4 disttocoast-ahead > 500 and item 5 disttocoast-ahead > 500) ; that is, if some pos ahead are close to land
;  let old-heading heading
;  facexy (item 0 (item 2 pos-list-daily)) (item 1 (item 2 pos-list-daily))
;  let tmp-heading heading
;  set heading old-heading
;  let tmp-angle subtract-headings tmp-heading old-heading  ; how much should I turn to go towards place I visited 2 days ago?
;  ifelse (tmp-angle < 0 ) [ rt -90 ] [ rt 90 ] ; turn left if that gets you back where you were previously


;  if (debug = 7 and (who = 0 or who = 1) ) [
;    print ""
;    print word "disttocoast-here: " [ disttocoast ] of patch-here
;    print word "disttocoast-ahead: " disttocoast-ahead
;  ]

  ; Stay on current dist from land if 1-4 km from land, or try to get there
  let dtc-diff 9999999 ; (list 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  let dtc-max -9999
  let dtc-min 9999
  let dtc-nbr-alike 9
  let dtc-nbr-max 9
  let dtc-nbr-min 9
  let tmp 99999
  let dtc-here [ disttocoast ] of patch-here

  let i 0
  while [ i < length ( disttocoast-ahead ) ] [
    set tmp abs( item i disttocoast-ahead  - dtc-here )
    if (tmp <= dtc-diff) [ ; if turning this much gets you to stay on same depth, choose that angle number
      set dtc-diff tmp
      set dtc-nbr-alike i
    ]
    ; set tmp item i disttocoast-ahead
    if (tmp <= dtc-min) [
      set dtc-min tmp
      set dtc-nbr-min i
    ]
    if (tmp > dtc-max) [
      set dtc-max tmp
      set dtc-nbr-max i
    ]

    set i i + 1
  ]

  let angles (list -80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80)
  let sel-angle -9999
  if ( [disttocoast] of patch-here > 4000 ) [ set sel-angle item dtc-nbr-min angles ]
  if ( [disttocoast] of patch-here < 1000 ) [ set sel-angle item dtc-nbr-max angles ]
  if ( [disttocoast] of patch-here <= 4000 and [disttocoast] of patch-here >= 2000) [ set sel-angle item dtc-nbr-alike angles ]
  set heading heading + sel-angle

;  if (debug = 7 and (who = 0 or who = 1) ) [
;    print word "dtc-here: " dtc-here
;    ; print disttocoast-ahead
;    ; print word "angles: " angles
;    print word "selected angle: " sel-angle
;    ; print ""
;  ]


  ; ### when to stop dispersing:  ###

  if ( not enough-water-ahead? ) [    ; distances are measured in number of cells here
    set disp-type 0
    if (debug = 7 and who = 0 ) [
      print "NO WATER, stop dispersing (2 - not enough water ahead) (porp 0)"
      pen-up
      set color red
    ]
    if (debug = 7 and who = 1 ) [
      print "NO WATER, stop dispersing (2 - not enough water ahead) (porp 1)"
      pen-up
      set color orange
    ]
  ]


  if ( not ( [ bathymetry ] of patch-ahead ( mean-disp-dist / 0.4 ) > min-depth ) ) [    ; distances are measured in number of cells here
    set disp-type 0
    if (debug = 7 and who = 0 ) [
      show "LOW water, stop dispersing (2) "
      pen-up
      set color red
    ]
    if (debug = 7 and who = 1 ) [
      show "LOW water, stop dispersing (2)"
      pen-up
      set color orange
    ]
  ]


  if ( disp-type = 2 ) [
    fd ( mean-disp-dist / 0.4 )
    ; set energy-level energy-level - 0.001 * e-use-per-km * mean-disp-dist / 0.4 ; REMOVED BY CARA
  ]
end ; end porp-disp-2


to porp-upd-mortality
  ; Introducing maximum age
  if (age > 30) [
    set list-of-dead-age lput (floor age) list-of-dead-age
    set list-of-dead-day lput (floor sim-day) list-of-dead-day
    die
  ]

  ; Mortality due to by-catch
  let daily-survival-prob exp( ln((1 - bycatch-prob)) / 360 )  ; Ok that only divided by 360, called once per day
  if (random-float 1 > daily-survival-prob) [
    set list-of-dead-age lput (floor age) list-of-dead-age
    set list-of-dead-day lput (floor sim-day) list-of-dead-day
    die
  ]

  ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;

  ; Prob of dying increases with decreasing body condition
 let m-mort-prob-const 10 ^ (2.176e-02 * x-surv-prob-const + -3.875e-04)                     ; EQN 63: obtained from fitting log-transformed possible values with a linear model in R

  let yearly-surv-prob (1 - (m-mort-prob-const * exp(- storage-level * x-surv-prob-const) )) ; EQN 62: Yearly survival probability
  if yearly-surv-prob < 0 [set yearly-surv-prob 0]

  let step-surv-prob 0
  ifelse (storage-level > 0.05 and yearly-surv-prob > 0 )                                    ; daily survival probability
  [ set step-surv-prob exp( ln(yearly-surv-prob) / (360) )]                                  ; EQN 64
  [ set step-surv-prob 0]

  if ( random-float 1 >= step-surv-prob ) [
    if (not with-lact-calf? or storage-level <= 0) [
      set list-of-dead-age lput (floor age) list-of-dead-age
      set list-of-dead-day lput (floor sim-day) list-of-dead-day
      if (debug = 10) [print word who " died of low body condition-Upd Mortality"]
     die
  ]]

  ; check for calf death
    if (with-lact-calf? = TRUE) [
      let storage-level-calf ((mass-calf - mass-struct-calf) / mass-calf)
      let yearly-surv-prob-calf (1 - (m-mort-prob-const * exp(- storage-level-calf * x-surv-prob-const) ))
      if yearly-surv-prob-calf < 0 [set yearly-surv-prob-calf 0]
      let step-surv-prob-calf 0
      ifelse (storage-level-calf > 0.05 and yearly-surv-prob-calf > 0) [ set step-surv-prob-calf exp( ln(yearly-surv-prob-calf) / 360 )] [set step-surv-prob-calf 0] ; if calf's storage level is above minimum storage level (0.05) calculate survival probability for the timestep
      if ( random-float 1 >= step-surv-prob-calf ) [
          set list-of-dead-age-calves lput (floor dsg-birth / 360) list-of-dead-age-calves
          set list-of-dead-day-calves lput (floor sim-day) list-of-dead-day-calves
          if (debug = 9) [print word who "'s calf died of low body condition-Upd Mortality"]
          ; reset mother lactation variables
          set m-BMR-calf 0
          set lgth-calf 0
          set mass-calf 0
          set blub-calf 0
          set mass-struct-calf 0
          set m-thermo-calf 0
          set max-grow-calf 0
          set m-growth-calf 0
          set m-blub-calf 0
          set e-calf 0
          set m-lact 0
          set m-lact-real 0
          set vital-costs 0
          set vital-costs-calf 0
          set wean-scale-fact 1
          set with-lact-calf? false
          set dsg-birth -99
          set IR-record-calf 0
          set n-calf-lost n-calf-lost + 1
        ]
    ]
  ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;
end ; end porp-upd-mortality


; Statistics and plots
to my-update-plots ; update histograms and other plots
  ; food level and number of porpoises
  set-current-plot "population"
  set-current-plot-pen "N x10"
  plot 10 * ( count porps )
  set-current-plot-pen "total food"
  plot (sum [ food-level ] of patches with [ food-level > 0 ])/ 2
  set-current-plot-pen "E x100"
  plot 1000 * mean [ storage-level ] of porps


  ; Histogram of porpoise storage-levels
  set-current-plot "porpoise-storage-level"
  let a [precision storage-level 2] of porps
  let c [precision ((mass-calf - mass-struct-calf) / mass-calf) 2] of porps with [with-lact-calf? = TRUE and sex-calf = "female"]
  set a sentence c a
  histogram a

   ; Age-class distribution in porpoise population
  set-current-plot "age-distribution"
  let aa [round age] of porps
  let cc [round (dsg-birth / 360)] of porps with [with-lact-calf? = TRUE and sex-calf = "female"]
  set aa sentence cc aa
  histogram aa

  ; Total blubber mass of porpoises of four classes (calves, juveniles, lactating mature, and non-lactating mature) ; filtered by weight to be comparable to Read 1990 and McLellan et al. 2002
  set-current-plot "blubber-mass"
  set-current-plot-pen "Lact"
  if year = sc-year [ask porps with  [with-lact-calf? = TRUE and weight >= 59.3 and weight <= 68.3] [plotxy (0.75 + random-float 0.1) (v-blub * dens-blub)]]
  set-current-plot-pen "Nonlact"
  if year = sc-year [ask porps with [with-lact-calf? = FALSE and age > age-of-maturity and weight >= 50.5 and weight <= 60.3][plotxy (0.55 + random-float 0.1) (v-blub * dens-blub)]]
  set-current-plot-pen "Juvs"
  if year = sc-year [ask porps with [age < age-of-maturity and lgth > 1.17 and weight >= 33.6 and weight <= 44][plotxy (0.35 + random-float 0.1) (v-blub * dens-blub)]]
  set-current-plot-pen "Calves"
  if year = sc-year [
   ask porps with  [with-lact-calf? = TRUE and mass-calf >= 23.3 and mass-calf <= 34.3][plotxy (0.15 + random-float 0.1) (v-blub-calf * dens-blub)]
   ask porps with  [lgth < 1.17 and weight >= 23.3 and weight <= 34.3][plotxy (0.15 + random-float 0.1) (v-blub * dens-blub)]
  ]

  ; Total body mass to length relationship
  set-current-plot "mass-length"
  if year = sc-year [
    ask porps [plotxy lgth weight]
    if any? porps with [with-lact-calf? = TRUE] [ask porps with [with-lact-calf? = TRUE] [plotxy lgth-calf mass-calf]]
  ]

  ;;; Seasonal measures:
  let dayofyear sim-day - ((year - 1) * 360)

  ; Average site-specific blubber depths
  set-current-plot "blubber-depths"
  set-current-plot-pen "axillary-D"
  if year = sc-year [plotxy dayofyear (mean [axillary-D] of porps)]
  set-current-plot-pen "axillary-L"
  if year = sc-year [plotxy dayofyear (mean [axillary-L] of porps)]
  set-current-plot-pen "axillary-V"
  if year = sc-year [plotxy dayofyear (mean [axillary-V] of porps)]
  set-current-plot-pen "CrIDF-D"
  if year = sc-year [plotxy dayofyear (mean [CrIDF-D] of porps)]
  set-current-plot-pen "CrIDF-L"
  if year = sc-year [plotxy dayofyear (mean [CrIDF-L] of porps)]
  set-current-plot-pen "CrIDF-V"
  if year = sc-year [plotxy dayofyear (mean [CrIDF-V] of porps)]
  set-current-plot-pen "CaIDF-D"
  if year = sc-year [plotxy dayofyear (mean [CaIDF-D] of porps)]
  set-current-plot-pen "CaIDF-L"
  if year = sc-year [plotxy dayofyear (mean [CaIDF-L] of porps)]
  set-current-plot-pen "CaIDF-V"
  if year = sc-year [plotxy dayofyear (mean [CaIDF-V] of porps)]


  ; Field metabolic rates
  set-current-plot "field-metabolic-rates"
  set-current-plot-pen "NonRepro"
  if year = sc-year [plotxy dayofyear (mean [(m-tot * 48) / 1000000] of porps with [m-lact-real = 0 and m-preg = 0])]
  set-current-plot-pen "Repro"
  if year = sc-year [plotxy dayofyear (mean [(m-tot * 48) / 1000000] of porps with [m-lact-real > 0 or m-preg > 0])]

  ; Average energy intake
  set-current-plot "energy-intake"
  set-current-plot-pen "Adults"
  if year = sc-year [plotxy dayofyear (mean [(daily-food / 1000000)] of porps with [m-lact-real = 0 and m-preg = 0])]
  set-current-plot-pen "Repro"
  if year = sc-year [plotxy dayofyear (mean [(daily-food / 1000000)] of porps with [m-lact-real > 0 or m-preg > 0])]

end


; File I/O

to file-setup
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on? true
  while [go-on?] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [
      set go-on? false
      file-open outfile
    ]
    if (repl-counter = 30) [
      set go-on? false
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  ;file-print (" id day year utm-x utm-y block age energy") ; header line (space-separated). Note that "bathy" is non-standard.
   file-print (" id year day utm-x utm-y dist-to-coast length weight age storage-level") ; header line (space-separated). Note that "bathy" is non-standard.
end


to file-setup-2
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 30) [
      set go-on false
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print (" day utm-x utm-y block") ; header line (space-separated). Note that "bathy" is non-standard
end


to file-setup-3   ; energy budget tracked porpoise debugging file
  let repl-counter 0
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 30) [
      set go-on false
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print (" day year age length weight mass-struct storage-level volume-blubber calf-mass calf-mass-struct pregnancy-status with-lactating-calf? swimming-speed daily-food IR-record total-metabolism maintenance-costs thermoregulation-costs locomotion-costs pregnancy-costs lactation-costs growth-costs") ; header line (space-separated).
end


to file-setup-4    ; energy budget tracked porpoise debugging file
  let repl-counter 0
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 30) [
      set go-on false
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print (" day year month blubber-depth-mean blubber-depth-sd FMR-mean FMR-sd energy-intake-mean energy-intake-sd") ; header line (space-separated).
end


to file-write-line
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
 ask porps [
    file-write who ; not followed by CR
     file-write year
     file-write round sim-day ; not followed by CR
    file-write round(xcor * 400 + xllcorner) ; 240 km wide and 400 km tall non-wrapped landscape divided into 400 x 400 m cells
    file-write round(ycor * 400 + yllcorner)
    file-write [ disttocoast ] of patch-here ; CHANGED BY CARA
    file-write lgth
    file-write weight
    file-write age
    file-write storage-level ; CHANGED BY CARA
    file-print ""    ; followed by CR   file-print ""    ; followed by CR
  ]
  ; add entries

  file-close
end


to file-write-line-2
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
  ask porp 0 [
    file-write round sim-day ; not followed by CR
    file-write round(xcor * 400 + xllcorner) ; 240 km wide and 400 km tall non-wrapped landscape divided into 400 x 400 m cells
    file-write round(ycor * 400 + yllcorner)
    file-write [ block ] of patch-here
;    file-write lgth
;    file-write weight
;    file-write age
;    file-write energy-level
    file-print ""    ; followed by CR
  ]
  file-close
end


to file-write-line-3
  ; Ask porpoise to write variables that can be used to debug energy budget procedures:
  file-open outfile ; append
  ; add entries
  ask porp 0 [
    file-write round sim-day
    file-write year
    file-write age
    file-write lgth
    file-write weight
    file-write mass-struct
    file-write storage-level ; CHANGED BY CARA
    file-write v-blub
    file-write mass-calf
    file-write mass-struct-calf
    file-write pregnancy-status
    file-write with-lact-calf?
    file-write swim-speed
    file-write daily-food
    file-write IR-record
    file-write m-tot
    file-write m-BMR
    file-write m-thermo
    file-write m-loco
    file-write m-preg
    file-write m-lact
    file-write m-growth
    file-print ""    ; followed by CR
  ]
  file-close
end


to file-write-line-4
  ; Ask porpoise to write variables that can be used to debug energy budget procedures:
  file-open outfile ; append
  ; add entries
    file-write round sim-day
    file-write year
    file-write month
    file-write mean [CaIDF-L] of porps
    file-write standard-deviation [CaIDF-L] of porps
    file-write mean [(m-tot * 48) / 1000000] of porps
    file-write standard-deviation [(m-tot * 48) / 1000000] of porps
    file-write mean [(daily-food / 1000000)] of porps
    file-write standard-deviation [(daily-food / 1000000)] of porps
    file-print ""    ; followed by CR
  file-close
end

; Main

to setup
  clear-all
  clear-turtles
  clear-output
  clear-drawing   ; clears pendown tracks
  clear-all-plots
  reset-ticks
  reset-timer

  ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
  ; ADDED BY CARA - for scenario outuputs
  set daily-pop-list []
  set fourty-km-list []
  set storage-level-deter-list []
  set mean-calf-mass-deter-list []
  set mean-juv-mass-deter-list []
  set mean-calf-mass-list []
  set mean-juv-mass-list []
  ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

  setup-global-parameters ; name changed by CARA
  landsc-setup
  porps-setup
  update-block-values

  if ( write-data != "off" and write-data != "one-porp" and write-data != "one-porp-energy-debug" and write-data != "all-porps-energy-debug") [
    file-setup
    print word "Writing to: " outfile
  ]
  if ( write-data = "one-porp-energy-debug") [
    file-setup-3
    print word "Writing to: " outfile
  ]

  if ( write-data = "all-porps-energy-debug") [
    file-setup-4
    print word "Writing to: " outfile
  ]

  reset-ticks

  ; Debugging code:
  if ( write-data = "one-porp") [
    file-setup-2
    print word "Writing porp 0 data to: " outfile
  ]
  if (debug = 1) [
    print ""
    print "Debug avoidance behaviour near land (debug=1):"
  ]
  if (debug = 2) [
    print ""
    print "Write turning angles before/after react. to land (debug=2):"
  ]
  if (debug = 3) [
    print ""
    print "Debugging CRW related turns (mod >=1) (debug=3):"
  ]
  if ( (debug = 4) and (model >= 2) ) [
    print ""
    print "Debugging attraction vector (mod >=2) (debug=4):"
  ]
  if ( (debug = 5) and (model >= 2) ) [
    print ""
    print "Debugging attraction + deterrence vectors (mod >=2)  (debug=6):"
  ]
  if ( debug = 6 ) [
    print ""
    print "Showing mean energy per 40 x 40 km block (debug=6)"
    print "day \tSE Anh (16) \tGB-N (33)  \tNV Born (42) \tFehm (52)"
  ]
  if (debug = 7) [
    print "Debugging dispersal behaviour"
    inspect porp 0
    inspect porp 1
  ]
  if (debug = 8) [ print "Debugging deterrence from wind farms and ships" ]
  if (debug = 9) [ print "Debugging reprod. and mortality" ]
  if (debug = 10) [ print "Debugging energy budget procedures and energy related mortality" ]
  if (debug = 11) [ print "Debugging energy intake procedures" ]
  if (debug = 12) [ print "Debugging storage procedures" ]
end


to daily-tasks
    ; Things to do daily or less often
    my-update-plots
    ;  landsc-display  ; turned off for performance

     if write-data = "daily" [ ; must write deployment pos before moving in movement models
     ask ( porps ) [ file-write-line ]
     ]

    ; write data to file:
   if write-data = "one-porp" [ ; one porpoise (porp 0), daily
      ask ( porp 0 ) [ file-write-line-2 ]
    ]

    if ( not (prev-year = year ) or ( time-step = 1 ) ) [
      if write-data = "yearly" [
        ask ( porps ) [ file-write-line ]
      ]
    ]

     if write-data = "one-porp-energy-debug" [ ; one porpoise (porp 0), daily
      ask ( porp 0 ) [ file-write-line-3 ]
    ]

     if write-data = "all-porps-energy-debug" [ ; averages for all porpoises, daily
      file-write-line-4
    ]

    ; Update daily average energy level and corresponding positions for porps and use it to start/stop dispersing
    ask ( porps ) [
       set age age + (1 / 360)

       ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
       ; keep track of average energy intake - averaged over a week as lots of variablity between timesteps
       if (remainder sim-day 7 ) = 0 [
       set daily-food ((sum food-intake-list) / 7) * IR-to-EA
       set food-intake-list []
       ]
       ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;

      let s-list-lgt length storage-level-daily                                   ; CHANGED BY CARA to reference storage level rather than energy level
      set storage-level-daily remove-item ( s-list-lgt - 1 ) storage-level-daily  ; CHANGED BY CARA
      let s-mean storage-level-sum / 48  ; 48 half-hour steps per day             ; CHANGED BY CARA
      set storage-level-daily fput (precision s-mean 3) storage-level-daily       ; CHANGED BY CARA
      set pos-list-daily fput ( list xcor ycor ) pos-list-daily
      set pos-list-daily remove-item ( length pos-list-daily - 1) pos-list-daily

      if ( model >= 3 ) [
        if ( disp-type = 0 ) [ ; i.e. not dispersing,
        if ( ( item 0 storage-level-daily ) < ( item 1 storage-level-daily ) and ( item 1 storage-level-daily ) < ( item 2 storage-level-daily )  and ( item 2 storage-level-daily ) < ( item 3 storage-level-daily ) ) [; CHANGED BY CARA
            set disp-type 1 ; decreasing energy for three days
            porp-disp-target-select
          ]
        ]
        if ( disp-type > 0 ) [
          ; storage level higher than any of the previous seven days, stop dispersing;
          if ( item 0 storage-level-daily ) > max (list item 1 storage-level-daily item 2 storage-level-daily item 3 storage-level-daily item 4 storage-level-daily item 5 storage-level-daily item 6 storage-level-daily item 7 storage-level-daily ) [; CHANGED BY CARA
            set disp-type 0
            if (debug = 7 and who = 0) [
              ask ( porp 0 ) [ pen-up ]
              set color red
              print  "Food found, stop disp., (porp 0)"
            ]
            if (debug = 7 and who = 1) [
              ask ( porp 1 ) [ pen-up ]
              set color red
              print  "Food found, stop disp., (porp 1)"
            ]
            set disp-target list 0 0
          ]
        ]
        if (disp-type = 2) [
          ; Energy level higher last week, turn towards place visited 3 days ago
          if ( mean( list(item 0 storage-level-daily) (item 1 storage-level-daily) (item 3 storage-level-daily)) < mean( list(item 6 storage-level-daily) (item 7 storage-level-daily) (item 8 storage-level-daily) (item 9 storage-level-daily) ) ) [ ; CHANGED BY CARA
            facexy (item 0 (item 2 pos-list-daily)) (item 1 (item 2 pos-list-daily))
            ;set heading ( heading + random 16 - 8 )
            ;rt 180
            fd 1
            if (debug = 7 and (who = 0 or who = 1)) [ show "RETURNING to prev pos" ]
          ]
        ]
      ] ; end model >=3 tasks

    ]


    ; Update parameters related to demography and blubber depth
    if ( model >= 4 ) [
      ask ( porps ) [
        porp-upd-mortality
        porp-upd-pregnancy-status

      ]

    ]

    ; if scenarios are running, collect outputs  ; ADDED bBY CARA FOR SCENARIOS
    if scenario != 0 [scenario-outputs]

    if ( debug = 6 ) [
      ; get mean energy of porps in blocks 16 (SE of Anholt), 33 (N Great Belt), 42 (NCV of Bornholm) and 52 (Fehmarn)
      type word round(sim-day) "\t"
      let s-list [ storage-level ] of porps with [ block = 16 ]        ; CHANGED BY CARA
      if ( length s-list > 0 ) [ type word (precision mean s-list 2) "      \t" ]
      if ( length s-list = 0 ) [ type word "- " "\t\t" ]
      set s-list [ storage-level ] of porps with [ block = 33 ]        ; CHANGED BY CARA
      if (length s-list > 0) [ type word (precision mean s-list 2) "      \t" ]
      if (length s-list = 0) [ type word "- " "\t\t" ]
      set s-list [ storage-level ] of porps with [ block = 42 ]        ; CHANGED BY CARA
      if (length s-list > 0) [ type word (precision mean s-list 2) "      \t" ]
      if (length s-list = 0) [ type word "- " "\t\t" ]
      set s-list [ storage-level ] of porps with [ block = 52 ]        ; CHANGED BY CARA
      if (length s-list > 0) [ print word (precision mean s-list 2) "      \t" ]
      if (length s-list = 0) [ print word "- " "\t\t" ]
    ]

end ; end daily-tasks

to calc-age-distrib
  ; age list at start of year
  set age-distrib [ ]
  let aa [round age] of porps
  let cc [round (dsg-birth / 360)] of porps with [with-lact-calf? = TRUE and sex-calf = "female"]   ; added to include female calves
  set aa sentence cc aa
  set age-distrib lput ( length filter [ ?1 -> ?1 = 0 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 1 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 2 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 3 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 4 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 5 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 6 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 7 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 8 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 9 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 10 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 11 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 12 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 13 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 14 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 15 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 16 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 17 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 18 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 19 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 20 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 21 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 22 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 23 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 = 24 ] aa ) age-distrib
  set age-distrib lput ( length filter [ ?1 -> ?1 >= 25 ] aa ) age-distrib

end ; end calc-age-distrib

to calc-mort-prob
  ; age distribution in list of dead
  let mort-age-distrib [ ]
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 0 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 1 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 2 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 3 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 4 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 5 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 6 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 7 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 8 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 9 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 10 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 11 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 12 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 13 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 14 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 15 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 16 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 17 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 18 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 19 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 20 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 21 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 22 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 23 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 = 24 ] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [ ?1 -> ?1 >= 25 ] list-of-dead-age ) mort-age-distrib

  ; report results
  if (debug = 0 and model >= 4) [
    if (sim-day < 2) [
      print " "
      print "Age class distrib. and yearly mort.:"
      print " y age n mort"
    ]
    foreach [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25] [ ?1 ->
      write round year - 1
      write ?1
      write (item ?1 age-distrib)
      print word " " (item ?1 mort-age-distrib)
    ]
  ]

  ; reset mort lists every year, after reporting:
  set list-of-dead-age [ ]
  set list-of-dead-day [ ]
end ; end calc-mort-prob

to yearly-tasks
  ask ( porps ) [ set mating-day round( random-normal (7.5 * 360 / 12) 20 )]
  ; Make list of number of dead per age class and corresp pop nos in prev year
  if (sim-day < 2) [ calc-age-distrib ]
  calc-mort-prob
  calc-age-distrib

end ; end yearly-tasks


to go
  if ( debug = "profile" ) [
    profiler:start
  ]
  if ( time-step = 0 ) [
    reset-timer
  ]

  if (model >= 3) [  ; deterrence behaviour
    porps-upd-deter
    if any? ss-ships [ships-deter-porps]
  ]

  ; for debugging:
  let EA 0

  ask porps [
    porp-move  ; This is the important step!
    if ( disp-type = 1 ) [ porp-disp-1 ] ; long-dist dispersal, set in daily-tasks
    if ( disp-type = 2 ) [ porp-disp-2 ] ; dispersal away from old pos, along coast
    if (debug = 8) [     ; if no longer scared by turbines, change color to std
    if (deter-strength = 0) [ ifelse(who = 0) [set color red] [set color orange]]
    ]
    ;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;
    ; this is where energy budget proceedures are called
    calc-swim-speed
    patch-zero-check
    ; assimilate energy
    energy-intake
    set EA e-assim
    ; then allocate it to:
    maintenance
    thermoregulation
    locomotion
    reproduction
    growth
    storage
    upd-state-variables
    ;;; ********************* ^^^ ADDED BY CARA ^^^ ********************** ;;;
]

  if ( model >= 3 and incl-ships?) [
    if scenario != 0 [
      if time-step = ((sc-year - 1)* 360 * 48) [ ships-setup ]; scenario setup-ADDED BY CARA FOR SCENARIOS
      ships-move ]]

  if (debug = 8) [ ; write porp xy near wind farms
    if (area = "Homogeneous" and wind-farms = "Line" and write-data = "off") [ file-noise-debug ]
  ]

  set time-step time-step + 1
  let prev-day floor sim-day
  set prev-month month
  set prev-quarter quarter
  set prev-year year
  set sim-day time-step / 48                                          ; day since start of sim
  set year floor ( sim-day / 360 ) + 1                                ; this gives 90 days per quarter
  set month ceiling ( (sim-day - (year - 1) * 360) / 30 )
  set quarter 1 + floor ( (sim-day + 30 - (year - 1) * 360) / 90 )    ; quarter 1 is Dec-Febr etc.
  if quarter = 5 [ set quarter 1 ]
  ;tick                                                               ; slows things down a lot

  ; Update food growth map and food level
  if ( not (prev-quarter = quarter )  ) [
    landsc-upd-maxent-level-map
  ]

  if ( remainder sim-day food-upd-interval ) = 0 [
    landsc-upd-food
  ]

  ; Things to do yearly:
  if ( not (prev-year = floor year ) or ( time-step = 1 ) ) [
    yearly-tasks
  ]

  ; Things to do daily (or less often):
  if ( not (prev-day = floor sim-day ) or ( time-step = 1 ) ) [
    daily-tasks
  ]

  if ( not (prev-month = month ) or ( time-step = 1 ) ) [
    if (write-data = "monthly") [
      ask ( porps ) [ file-write-line ]
    ]
    temp-sal-update     ; update temperature and salinity maps monthly ; ADDED BY CARA
  ]

  ; Stop running and print summary
  if ( time-step / 48 = max-sim-day ) [
    let ttt word "Time (" time-step
    set ttt word ttt " half-hour intervals): "
    set ttt word ttt timer
    print word ttt " sec"
    ; print ( list "Time (" max-tick "ticks): " timer "sec" )
    if (write-data != "off") [ file-close  ]
    stop
  ]

  if ( debug = "profile" ) [
    profiler:stop
    print profiler:report
    profiler:reset
  ]

   if (debug = 10) [ if (remainder sim-day 7 ) = 0 [
    print "Debugging energy budget calcs: all animals - weekly"
    ; maintenance
    print (word "Maintenance       " "min: " (precision (min [m-BMR] of porps) 1) " mean: " (precision (mean [m-BMR] of porps) 1) " max: " (precision (max [m-BMR] of porps) 1))
    ; thermoregulation
    print (word "Thermal costs     " "min: " (precision (min [m-thermo] of porps) 1) " mean: " (precision (mean [m-thermo] of porps) 1) " max: " (precision (max [m-thermo] of porps) 1))
    ; locomotion
    print (word "Locomotion        " "min: " (precision (min [m-loco] of porps) 1) " mean: " (precision (mean [m-loco] of porps) 1) " max: " (precision (max [m-loco] of porps) 1))
    ; pregnancy
    if any? porps with [m-preg > 0] [print (word "Pregnancy         " "min: " (precision (min [m-preg] of porps with [m-preg > 0]) 1) " mean: " (precision (mean [m-preg] of porps with [m-preg > 0]) 1) " max: " (precision (max [m-preg] of porps with [m-preg > 0]) 1))]
    ;lactation
    if any? porps with [m-lact > 0] [print (word "Lactation         " "min: " (precision (min [m-lact] of porps with [m-lact > 0]) 1) " mean: " (precision (mean [m-lact] of porps with [m-preg > 0]) 1) " max: " (precision (max [m-lact] of porps with [m-lact > 0]) 1))]
    ; growth
    print (word "Growth            " "min: " (precision (min [m-growth] of porps) 1) " mean: " (precision (mean [m-growth] of porps) 1) " max: " (precision (max [m-growth] of porps) 1))
    ; calf costs
    if any? porps with [m-lact > 0] [print (word "Calf maintenance " "min: " (precision (min [ m-BMR-calf] of porps with [m-lact > 0]) 1) " mean: " (precision (mean [ m-BMR-calf] of porps with [m-preg > 0]) 1) " max: " (precision (max [ m-BMR-calf] of porps with [m-lact > 0]) 1))]
    if any? porps with [m-lact > 0] [print (word "Calf thermal     " "min: " (precision (min [m-thermo-calf] of porps with [m-lact > 0]) 1) " mean: " (precision (mean [m-thermo-calf] of porps with [m-preg > 0]) 1) " max: " (precision (max [m-thermo-calf] of porps with [m-lact > 0]) 1))]
    if any? porps with [m-lact > 0] [print (word "Calf growth      " "min: " (precision (min [m-growth-calf] of porps with [m-lact > 0]) 1) " mean: " (precision (mean [m-growth-calf] of porps with [m-preg > 0]) 1) " max: " (precision (max [m-growth-calf] of porps with [m-lact > 0]) 1))]
    if any? porps with [m-lact > 0] [print (word "calf blubber     " "min: " (precision (min [m-blub-calf] of porps with [m-lact > 0]) 1) " mean: " (precision (mean [m-blub-calf] of porps with [m-preg > 0]) 1) " max: " (precision (max [m-blub-calf] of porps with [m-lact > 0]) 1))]
    print "END: Debugging energy budget calcs: all animals"
    print ""
  ]]

  if (debug = 11) [ if (remainder sim-day 7 ) = 0 [
   print "Debugging energy intake: all animals - weekly"
   print (word "IR-Record range    " "min: " (precision (min [IR-Record] of porps) 3) " mean: " (precision (mean [IR-Record] of porps) 3) " max: " (precision (max [IR-Record] of porps) 3))
   print (word "daily-food range   " "min: " (precision (min [daily-food / IR-to-EA] of porps) 3) " mean: " (precision (mean [daily-food / IR-to-EA] of porps) 3) " max: " (precision (max [daily-food / IR-to-EA] of porps) 3))
   print "END: Debugging energy intake: all animals"
   print ""
  ]]


    if (debug = 12) [ if (remainder sim-day 7 ) = 0 [
    print "Debugging storage calcs: all animals - weekly"
    ; storage levels
    print (word "Storage levels         " "min: " (precision (min [storage-level] of porps) 3) " mean: " (precision (mean [storage-level] of porps) 3) " max: " (precision (max [storage-level] of porps) 3))
    ; blubber volume
    print (word "Blubber volume         " "min: " (round (min [v-blub] of porps)) " mean: " (round (mean [v-blub] of porps)) " max: " (round (max [v-blub] of porps)))
    ; calf storage levels
    if any? porps with [m-lact > 0] [print (word "Calf storage levels    " "min: " (precision (min [((mass-calf - mass-struct-calf) / mass-calf)] of porps) 3) " mean: " (precision (mean [((mass-calf - mass-struct-calf) / mass-calf)] of porps) 3) " max: " (precision (max [((mass-calf - mass-struct-calf) / mass-calf)] of porps) 3))]
    ; calf blubber volume
    if any? porps with [m-lact > 0] [print (word "Calf blubber volume    " "min: " (round (min [v-blub-calf] of porps)) " mean: " (round (mean [v-blub-calf] of porps)) " max: " (round (max [v-blub-calf] of porps)))]

    ; site specific blubber depths
    print (word "Blubber depth: AxD     " "min: " (precision (min [axillary-D] of porps) 2) " mean: " (precision (mean [axillary-D] of porps) 2) " max: " (precision (max [axillary-D] of porps) 2))
    print (word "Blubber depth: AxL     " "min: " (precision (min [axillary-L] of porps) 2) " mean: " (precision (mean [axillary-L] of porps) 2) " max: " (precision (max [axillary-L] of porps) 2))
    print (word "Blubber depth: AxV     " "min: " (precision (min [axillary-V] of porps) 2) " mean: " (precision (mean [axillary-V] of porps) 2) " max: " (precision (max [axillary-V] of porps) 2))
    print (word "Blubber depth: CrIDF-D " "min: " (precision (min [CrIDF-D] of porps) 2) " mean: " (precision (mean [CrIDF-D] of porps) 2) " max: " (precision (max [CrIDF-D] of porps) 2))
    print (word "Blubber depth: CrIDF-L " "min: " (precision (min [CrIDF-L] of porps) 2) " mean: " (precision (mean [CrIDF-L] of porps) 2) " max: " (precision (max [CrIDF-L] of porps) 2))
    print (word "Blubber depth: CrIDF-V " "min: " (precision (min [CrIDF-V] of porps) 2) " mean: " (precision (mean [CrIDF-V] of porps) 2) " max: " (precision (max [CrIDF-V] of porps) 2))
    print (word "Blubber depth: CaIDF-D " "min: " (precision (min [CaIDF-D] of porps) 2) " mean: " (precision (mean [CaIDF-D] of porps) 2) " max: " (precision (max [CaIDF-D] of porps) 2))
    print (word "Blubber depth: CaIDF-L " "min: " (precision (min [CaIDF-L] of porps) 2) " mean: " (precision (mean [CaIDF-L] of porps) 2) " max: " (precision (max [CaIDF-L] of porps) 2))
    print (word "Blubber depth: CaIDF-V " "min: " (precision (min [CaIDF-V] of porps) 2) " mean: " (precision (mean [CaIDF-V] of porps) 2) " max: " (precision (max [CaIDF-V] of porps) 2))
    print "END: Debugging storage calcs: all animals"
    print ""
  ]]

end


;;; ********************* vvv ADDED BY CARA vvv ********************** ;;;



to setup-global-parameters  ; setup global parameters

  ;;; Energy budget parameters
  ; Submodel: Energy Intake
   set IR-coef 0.0004              ; Ingestion rate coefficient - estimated by ensuring IR-record sizes are not too large
   set satiation-c 10              ; Satiation constant - estimated using population storage levels
   set AE-food  0.82               ; Assimilation efficiency of food - Kriete 1995 -  killer whales fed fish diet
   set IR-to-EA 113750000          ; Ingested food to energy available - calibrated

  ; Submodel: LOCOMOTION
   set prop-eff 0.81               ; Propeller efficiency - Fish 1993

  ; Submodel: REPRODUCTION
   set age-of-maturity 3.44        ; Age of maturity - Read (1990) and Caswell
   set max-mass-f 8                ; Max mass of fetus - Lockyer and Kinze 2003
   set f-growth-c 0.0066858          ; Fetal growth constant - calculated from Lockyer & Kinze 2003 relationship
   set percent-lip-f 0.285         ; Fetal blubber percent composition - Blubber percent from McLellan et al. 2002 and percent lipid of blubber as 68.3% for neonates from Lockyer 1995 ( Marine Mammales: Biology and Conservation pg 111)
   set percent-pro-f 0.139         ; Fetal protein percent composition - LM percent from McLellan et al. 2002 and percent protein of LM from Lockyer 1991 (sperm whale)
   set lact-eff 0.84               ; Efficiency of producing milk - Anderson and Fedak 1987 (grey seal)
   set repro-min-SL 0.10           ; Minimum storage level for reproductive energy allocation - Beltran et al. 2017
   set calf-idl-SL 0.375           ; Mean percent blubber for calves in McLellan et al. 2002
   set t-gest 300                  ; Gestation period - Lockyer et al., 2003
   set t-nurs 240                  ; Lactation period - Lockyer et al., 2003; Lockyer and Kinze, 2003
   set pregnancy-rate 0.67         ; Pregnancy rate - Sørensen & Kinze 1994

  ; Submodel: GROWTH
   set ED-lip 39.5 * 1000000       ; Lipid energy density - Brody 1968, Blaxter 1989, Worthy 1982
   set ED-pro 23.6 * 1000000       ; Protein energy density - Brody 1968

  ; Submodel: STORAGE
   set dens-blub 0.00092           ; Density of lipid - Parry 1949 adjusted for water content

  ; Submodel: LIFE HISTORY
   set x-surv-prob-const 13.5      ; Survival probability constanct - calibrated

  ;;; Environmental parameters - Slider params from Nabe-Nielsen et al. 2014 not changed in this version of the model
  ; dispersal params:
  set min-disp-depth 4.0           ; in m
  set n-disp-targets 12
  set mean-disp-dist 1.6           ; km / 30 min
  set min-dist-to-target 100       ; km

  set maxU 1.00
  set food-growth-rate 0.10        ; rU
  set gravity 9.8

end


to porps-setup-params

 ; Submodel: MAINTENANCE
  set B0 11.13 + random-float 0.04 - random-float 0.04               ; Maintenance normalization constant - calibrated

  ; Submodel: THERMOREGULATION
  set for-convec-scal-coef-list []                                   ; List setup to hold cone forced convection scaling coefficient values
  set hC-low-lim-coef-list []                                        ; List setup to hold cone heat transfer lower limit values
  set T-core random-normal 36.7 0.4                                  ; Core temperature - Blanchet, Wahlberg, and Acquarone 2008
  set kB random-normal 0.10 0.01                                     ; Thermal conductivity of blood free blubber - Worthy and Edwards 1990

  ; Submodel: LOCOMOTION
   set lambda random-normal 0.25 0.016                               ; Ratio of active to passive drag - calibrated

  ; Submodel: REPRODUCTION
  set m-preg 0                                                       ; Initialize pregnancy and lactation costs
  set m-lact 0
  set preg-chance 0                                                  ; Initialize pregnancy chance as 0

  ; Submodel: GROWTH
  set struct-mass-perc-pro random-normal 0.2629 0.0086               ; Percent structural mass that is protein - Lockyer 1991 (sperm whale)
  set struct-mass-perc-lip random-normal 0.0288 0.0114               ; Percent structural mass that is lipid - Lockyer 1991 (sperm whale)
  set DE-lip 0.74 + random-float 0.16                                ; Deposition efficiency of lipid - Malavear 2002, Pullar & Webster 1977
  set DE-pro 0.43 + random-float 0.13                                ; Deposition efficiency of protein - Malavear 2002, Pullar & Webster 1977
  set ED-lean-mass ((struct-mass-perc-pro * ED-pro ) + (struct-mass-perc-lip * ED-lip))                                                 ; EQN 55: Lean mass energy density
  set DE-lean-mass ((struct-mass-perc-pro * DE-pro ) + (struct-mass-perc-lip * DE-lip))/ (struct-mass-perc-pro + struct-mass-perc-lip)  ; EQN 56: Lean mass deposition efficiency
  set m-str-k random-normal 1.16 0.33                                ; Mass growth constant for females from Galatius & Kinze unps. data fitted with a VB curve
  set m-str-inf random-normal 46.68 3.86                             ; Mass asymptotic value for females from Galatius & Kinze unps. data fitted with a VB curve
  set lgth-inf random-normal 158.12 4.68                             ; Length asymptotic value for females from Galatius & Kinze unps. data fitted with a gompertz curve
  set lgth-0 random-normal 94.82 1.69                                ; Length initial value for females from Galatius & Kinze unps. data fitted with a gompertz curve
  set lgth-k random-normal 0.41 0.06                                 ; Length k value from Galatius & Kinze unps. data fitted with a gompertz curve

  ; Submodel: STORAGE
  set perc-lip-blub random-normal 0.816 0.036                        ; Percent blubber that is lipid - Worthy and Edwards 1990


; initialize morphometrics based on age - For all groups: 0, 1, 2, 3, 4, 5, 6, 7 & up
 ifelse age = 0
   [ let porp-init-0 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Zero.csv"  ; Import csv
     let porp-init-matrix-0 matrix:from-row-list porp-init-0                                               ; Create matrix from csv
     let n-porps-0 matrix:get-column porp-init-matrix-0 0                                                  ; Pull the first column as a list to get max row number
     let n-0 (1 + (random (length n-porps-0 - 1)))                                                         ; Generate random row number
     set lgth matrix:get porp-init-matrix-0 n-0 1                                                          ; Pull length from csv
     set weight matrix:get porp-init-matrix-0 n-0 2                                                        ; Pull mass from csv
     set v-blub (matrix:get porp-init-matrix-0 n-0 3) / dens-blub                                          ; Set volume of blubber as the blubber mass multiplied by the blubber density in cm3
     set mass-struct (matrix:get porp-init-matrix-0 n-0 2) - (matrix:get porp-init-matrix-0 n-0 3)         ; Set structural mass as the difference between weight and blubber mass in kg

     set axillary-D matrix:get porp-init-matrix-0 n-0 8                                                    ; Axillary dorsal site blubber depth in cm
     set axillary-L matrix:get porp-init-matrix-0 n-0 9                                                    ; Axillary lateral site blubber depth in cm
     set axillary-V matrix:get porp-init-matrix-0 n-0 10                                                   ; Axillary ventral site blubber depth in cm
     set CrIDF-D matrix:get porp-init-matrix-0 n-0 11                                                      ; Cranial insertion of the dorsal fin dorsal site blubber depth in cm
     set CrIDF-L matrix:get porp-init-matrix-0 n-0 12                                                      ; Cranial insertion of the dorsal fin lateral site blubber depth in cm
     set CrIDF-V matrix:get porp-init-matrix-0 n-0 13                                                      ; Cranial insertion of the dorsal fin ventral site blubber depth in cm
     set CaIDF-D matrix:get porp-init-matrix-0 n-0 14                                                      ; Caudal insertion of the dorsal fin dorsal site blubber depth in cm
     set CaIDF-L matrix:get porp-init-matrix-0 n-0 15                                                      ; Caudal insertion of the dorsal fin lateral site blubber depth in cm
     set CaIDF-V matrix:get porp-init-matrix-0 n-0 16 ]                                                    ; Caudal insertion of the dorsal fin ventral site blubber depth in cm

   [ ifelse age = 1                                                                                        ; 1 year old initilization
     [ let porp-init-1 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_One.csv"
       let porp-init-matrix-1 matrix:from-row-list porp-init-1
       let n-porps-1 matrix:get-column porp-init-matrix-1 0
       let n-1 (1 + (random (length n-porps-1 - 1)))
       set lgth matrix:get porp-init-matrix-1 n-1 1
       set weight matrix:get porp-init-matrix-1 n-1 2
       set v-blub (matrix:get porp-init-matrix-1 n-1 3) / dens-blub
       set mass-struct (matrix:get porp-init-matrix-1 n-1 2) - (matrix:get porp-init-matrix-1 n-1 3)

       set axillary-D matrix:get porp-init-matrix-1 n-1 8
       set axillary-L matrix:get porp-init-matrix-1 n-1 9
       set axillary-V matrix:get porp-init-matrix-1 n-1 10
       set CrIDF-D matrix:get porp-init-matrix-1 n-1 11
       set CrIDF-L matrix:get porp-init-matrix-1 n-1 12
       set CrIDF-V matrix:get porp-init-matrix-1 n-1 13
       set CaIDF-D matrix:get porp-init-matrix-1 n-1 14
       set CaIDF-L matrix:get porp-init-matrix-1 n-1 15
       set CaIDF-V matrix:get porp-init-matrix-1 n-1 16 ]

       [ ifelse age = 2
          [ let porp-init-2 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Two.csv"
            let porp-init-matrix-2 matrix:from-row-list porp-init-2
            let n-porps-2 matrix:get-column porp-init-matrix-2 0
            let n-2 (1 + (random (length n-porps-2 - 1)))
            set lgth matrix:get porp-init-matrix-2 n-2 1
            set weight matrix:get porp-init-matrix-2 n-2 2
            set v-blub (matrix:get porp-init-matrix-2 n-2 3) / dens-blub
            set mass-struct (matrix:get porp-init-matrix-2 n-2 2) - (matrix:get porp-init-matrix-2 n-2 3)

            set axillary-D matrix:get porp-init-matrix-2 n-2 8
            set axillary-L matrix:get porp-init-matrix-2 n-2 9
            set axillary-V matrix:get porp-init-matrix-2 n-2 10
            set CrIDF-D matrix:get porp-init-matrix-2 n-2 11
            set CrIDF-L matrix:get porp-init-matrix-2 n-2 12
            set CrIDF-V matrix:get porp-init-matrix-2 n-2 13
            set CaIDF-D matrix:get porp-init-matrix-2 n-2 14
            set CaIDF-L matrix:get porp-init-matrix-2 n-2 15
            set CaIDF-V matrix:get porp-init-matrix-2 n-2 16 ]

            [ ifelse age = 3
              [ let porp-init-3 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Three.csv"
                let porp-init-matrix-3 matrix:from-row-list porp-init-3
                let n-porps-3 matrix:get-column porp-init-matrix-3 0
                let n-3 (1 + (random (length n-porps-3 - 1)))
                set lgth matrix:get porp-init-matrix-3 n-3 1
                set weight matrix:get porp-init-matrix-3 n-3 2
                set v-blub (matrix:get porp-init-matrix-3 n-3 3) / dens-blub
                set mass-struct (matrix:get porp-init-matrix-3 n-3 2) - (matrix:get porp-init-matrix-3 n-3 3)

                set axillary-D matrix:get porp-init-matrix-3 n-3 8
                set axillary-L matrix:get porp-init-matrix-3 n-3 9
                set axillary-V matrix:get porp-init-matrix-3 n-3 10
                set CrIDF-D matrix:get porp-init-matrix-3 n-3 11
                set CrIDF-L matrix:get porp-init-matrix-3 n-3 12
                set CrIDF-V matrix:get porp-init-matrix-3 n-3 13
                set CaIDF-D matrix:get porp-init-matrix-3 n-3 14
                set CaIDF-L matrix:get porp-init-matrix-3 n-3 15
                set CaIDF-V matrix:get porp-init-matrix-3 n-3 16 ]

                [ ifelse age = 4
                  [ let porp-init-4 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Four.csv"
                    let porp-init-matrix-4 matrix:from-row-list porp-init-4
                    let n-porps-4 matrix:get-column porp-init-matrix-4 0
                    let n-4 (1 + (random (length n-porps-4 - 1)))
                    set lgth matrix:get porp-init-matrix-4 n-4 1
                    set weight matrix:get porp-init-matrix-4 n-4 2
                    set v-blub (matrix:get porp-init-matrix-4 n-4 3) / dens-blub
                    set mass-struct (matrix:get porp-init-matrix-4 n-4 2) - (matrix:get porp-init-matrix-4 n-4 3)

                    set axillary-D matrix:get porp-init-matrix-4 n-4 8
                    set axillary-L matrix:get porp-init-matrix-4 n-4 9
                    set axillary-V matrix:get porp-init-matrix-4 n-4 10
                    set CrIDF-D matrix:get porp-init-matrix-4 n-4 11
                    set CrIDF-L matrix:get porp-init-matrix-4 n-4 12
                    set CrIDF-V matrix:get porp-init-matrix-4 n-4 13
                    set CaIDF-D matrix:get porp-init-matrix-4 n-4 14
                    set CaIDF-L matrix:get porp-init-matrix-4 n-4 15
                    set CaIDF-V matrix:get porp-init-matrix-4 n-4 16 ]

                    [ ifelse age = 5
                      [ let porp-init-5 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Five.csv"
                        let porp-init-matrix-5 matrix:from-row-list porp-init-5
                        let n-porps-5 matrix:get-column porp-init-matrix-5 0
                        let n-5 (1 + (random (length n-porps-5 - 1)))
                        set lgth matrix:get porp-init-matrix-5 n-5 1
                        set weight matrix:get porp-init-matrix-5 n-5 2
                        set v-blub (matrix:get porp-init-matrix-5 n-5 3) / dens-blub
                        set mass-struct (matrix:get porp-init-matrix-5 n-5 2) - (matrix:get porp-init-matrix-5 n-5 3)

                        set axillary-D matrix:get porp-init-matrix-5 n-5 8
                        set axillary-L matrix:get porp-init-matrix-5 n-5 9
                        set axillary-V matrix:get porp-init-matrix-5 n-5 10
                        set CrIDF-D matrix:get porp-init-matrix-5 n-5 11
                        set CrIDF-L matrix:get porp-init-matrix-5 n-5 12
                        set CrIDF-V matrix:get porp-init-matrix-5 n-5 13
                        set CaIDF-D matrix:get porp-init-matrix-5 n-5 14
                        set CaIDF-L matrix:get porp-init-matrix-5 n-5 15
                        set CaIDF-V matrix:get porp-init-matrix-5 n-5 16 ]

                        [ ifelse age = 6
                           [ let porp-init-6 csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Six.csv"
                             let porp-init-matrix-6 matrix:from-row-list porp-init-6
                             let n-porps-6 matrix:get-column porp-init-matrix-6 0
                             let n-6 (1 + (random (length n-porps-6 - 1)))
                             set lgth matrix:get porp-init-matrix-6 n-6 1
                             set weight matrix:get porp-init-matrix-6 n-6 2
                             set v-blub (matrix:get porp-init-matrix-6 n-6 3) / dens-blub
                             set mass-struct (matrix:get porp-init-matrix-6 n-6 2) - (matrix:get porp-init-matrix-6 n-6 3)

                             set axillary-D matrix:get porp-init-matrix-6 n-6 8
                             set axillary-L matrix:get porp-init-matrix-6 n-6 9
                             set axillary-V matrix:get porp-init-matrix-6 n-6 10
                             set CrIDF-D matrix:get porp-init-matrix-6 n-6 11
                             set CrIDF-L matrix:get porp-init-matrix-6 n-6 12
                             set CrIDF-V matrix:get porp-init-matrix-6 n-6 13
                             set CaIDF-D matrix:get porp-init-matrix-6 n-6 14
                             set CaIDF-L matrix:get porp-init-matrix-6 n-6 15
                             set CaIDF-V matrix:get porp-init-matrix-6 n-6 16 ]


                           [ let porp-init-7u csv:from-file "porpoise-initialization-files/PorpoiseInitializationBlubber_Seven_up.csv"
                             let porp-init-matrix-7u matrix:from-row-list porp-init-7u
                             let n-porps-7u matrix:get-column porp-init-matrix-7u 0
                             let n-7u (1 + (random (length n-porps-7u - 1)))
                             set lgth matrix:get porp-init-matrix-7u n-7u 1
                             set weight matrix:get porp-init-matrix-7u n-7u 2
                             set v-blub (matrix:get porp-init-matrix-7u n-7u 3) / dens-blub
                             set mass-struct (matrix:get porp-init-matrix-7u n-7u 2) - (matrix:get porp-init-matrix-7u n-7u 3)

                             set axillary-D matrix:get porp-init-matrix-7u n-7u 8
                             set axillary-L matrix:get porp-init-matrix-7u n-7u 9
                             set axillary-V matrix:get porp-init-matrix-7u n-7u 10
                             set CrIDF-D matrix:get porp-init-matrix-7u n-7u 11
                             set CrIDF-L matrix:get porp-init-matrix-7u n-7u 12
                             set CrIDF-V matrix:get porp-init-matrix-7u n-7u 13
                             set CaIDF-D matrix:get porp-init-matrix-7u n-7u 14
                             set CaIDF-L matrix:get porp-init-matrix-7u n-7u 15
                             set CaIDF-V matrix:get porp-init-matrix-7u n-7u 16 ]

     ]]]]]]

  set age age + 0.5                                                                                                                  ; Increase age by 0.5 to account for Jan 1st start date
  ifelse pregnancy-status != 1 [set surface-area 0.093 * weight ^ 0.57] [set surface-area 0.093 * (weight + (2 * mass-f)) ^ 0.57]    ; Calculate surface area - Worthy and Edwards 1990
  set v-blub-calf 0                                                                                                                  ; Initialize calf blubber volume as 0
  set v-blub-min (weight * 0.05) / dens-blub                                                                                         ; Minimum blubber as 5% of weight
  set v-blub-Repro (weight * repro-min-SL / dens-blub)
  set e-repo-min (v-blub-Repro - v-blub-min) * dens-blub * ED-lip                                                                    ; Calculate energy required for reproductive threshold
  set SL-mean ((lgth * -0.3059) + 0.7066)* IR-temp-mod                                                                               ; Mean storage level percentages from McLellan et al. 2002 for mature, calf, and immature porpoises fit with a linear relationship to length and adjusted using temperature modifier
  set v-blub-mean (weight * SL-mean / dens-blub)                                                                                     ; Converted to blubber volume
  set v-blub-calf-idl (mass-calf * 0.375 / dens-blub)                                                                                   ; Mean percent blubber for calves in McLellan et al. 2002
  set e-storage (v-blub - v-blub-min) * dens-blub * ED-lip                                                                           ; Storage energy as the amount of energy stored in blubber over minimum threshold
  set storage-level (weight - mass-struct) / weight                                                                                  ; Storage level calculation

  upd-blubber-depths

  ; initialize tracked values as zero
  set daily-food 0
  set food-intake-list []
  set storage-level-sum 0
  set storage-level-daily ( list 0 0 0 0 0 0 0 0 0 0 )
  set IR-record 0
  set IR-record-calf 0
end




to temp-sal-update
    ; Update temperature and salinity maps

    if ( area = "Kattegat" ) [
    let mth month
    if mth = 0 [set mth 1]

    let yyy word "TemperatureData/TempData_Month" mth
    set yyy word yyy ".asc"
    set temp-data gis:load-dataset word path yyy
    gis:apply-raster temp-data temp-w
    ask patches [ ifelse not ( temp-w >= 0 ) [ set temp-w -1 ] [ set temp-w temp-w ] ]

    set yyy word "SalinityData/SalData_Month" mth
    set yyy word yyy ".asc"
    set salinity-data gis:load-dataset word path yyy
    gis:apply-raster salinity-data salinity-w
    ask patches [ ifelse not ( salinity-w >= 0 ) [ set salinity-w -1 ] [ set salinity-w salinity-w ] ]
    ]

  ; landsc-display

  ask water-patches [thermo-prop-update] ; Update thermophysical properties of seawater

  ; Update ingestion rate temperature modifier
  set mean-temp mean [temp-w] of water-patches
  let xxxxx ((1 / 6.366) * mean-temp) * (180 / pi)
  let yyyyy (cos xxxxx)
  set IR-temp-mod ((1 / 5)* yyyyy + 1)

end


to thermo-prop-update

; -------- Lookup tables for assigning thermophysical properties of water ------- ;

    ; setup density
    let coli position (round temp-w) row-dens
    let rowi position (round salinity-w) col-dens
    set density-w matrix:get density-matrix rowi coli

    ; setup thermal conductivity of water
    let coli1 position (round temp-w) row-conduct
    let rowi1 position (round salinity-w) col-conduct
    set conductivity-w matrix:get conduct-matrix rowi1 coli1

    ; setup coefficient of thermal expansion of water
    let coli2 position (round temp-w) row-ceote
    let rowi2 position (round salinity-w) col-ceote
    set coef-thermal-expansion-w matrix:get ceote-matrix rowi2 coli2

    ; setup specific heat of water
    let coli3 position (round temp-w) row-sp-heat
    let rowi3 position (round salinity-w) col-sp-heat
    let spec-heat-w matrix:get sp-heat-matrix rowi3 coli3

    ; setup dynamic viscocity
    let coli4 position (round temp-w) row-dyn-vis
    let rowi4 position (round salinity-w) col-dyn-vis
    let dynamic-vis-w matrix:get dynamic-vis-matrix rowi4 coli4

; -------- END: Lookup tables for assigning thermophysical properties of water ------- ;

    ; calculate prandtl number
    set Pr-w (dynamic-vis-w * spec-heat-w / conductivity-w)

    ; calculate the kinematic viscocity of water
    set kin-visc-w dynamic-vis-w / density-w
end


to patch-zero-check                                                                      ; some patches where porpoises move will have zero values. These should be filled the first time they are encountered.
  if [temp-w] of patch-here <= 0 [
  ask patch-here [
      set temp-w mean [temp-w] of water-patches in-radius 10
      set salinity-w mean [salinity-w] of water-patches in-radius 10
      thermo-prop-update
  ]
    set water-patches patches with [temp-w > 0 and salinity-w > 0 and bathymetry > 0]
]
  if [kin-visc-w] of patch-here = 0 and [temp-w] of patch-here > 0 [ ask patch-here [ thermo-prop-update ]]
end

;    --------------------------------------------
;    |      BEGIN ENERGY BUDGET PROCEDURES      |
;    --------------------------------------------


to energy-intake
  let IR-struct-mass ((IR-coef * (mass-struct ^ 0.75))  * IR-temp-mod)         ; EQN 13: Calculate max ingestion rate for tick based on mass
  let preg-IR-sup (m-preg / (AE-food * IR-to-EA))                              ; EQN 15: Pregnant females will increase their food intake to cover pregnancy costs - Rojano-Doñate et al. 2018
  let lact-IR-sup (m-lact-real /(AE-food * IR-to-EA))                          ; EQN 16: Lactating females will increase their food intake to cover lactation costs - Douhard et al. 2016

  ; Calculate the amount of food needed to be ingested by the calf based on its size - EQN 23
  let IR-struct-mass-calf 0
  if with-lact-calf? = true and wean-scale-fact < 1 [set IR-struct-mass-calf (((IR-coef * mass-struct-calf ^ 0.75) * (1.00 - wean-scale-fact)) * IR-temp-mod)]

  ; Add total ingestion rate for the timestep
  let IR-timestep IR-struct-mass + preg-IR-sup + lact-IR-sup                   ; EQN 17

  ; porpoises keep track of intake rate when they haven't encountered food to compensate - EQN 18
  set IR-record IR-record + IR-timestep
  if with-lact-calf? = TRUE and wean-scale-fact < 1 [ set IR-record-calf IR-record-calf + IR-struct-mass-calf ]

  ; Check food levels of patch here and adjust IR
  let IR-real 0
  let IR-real-calf  0
  let split 0
  let food-available [food-level] of patch (item 0 (item 1 pos-list)) (item 1 (item 1 pos-list))  ; item 0 pos-list was the last added element, i.e. the current position

  if (debug = 11) [ if food-available < 0 [ print (word "WARNING: Food-level of " (patch (item 0 (item 1 pos-list)) (item 1 (item 1 pos-list))) " < 0")]]

  ifelse with-lact-calf? = FALSE or wean-scale-fact = 1
  [ ifelse food-available > IR-record
    [set IR-real IR-record]
    [set IR-real food-available]
  ]
  [ ifelse food-available > IR-record + IR-record-calf
    [
      set IR-real IR-record
      set IR-real-calf IR-record-calf
    ]
    [ set split ( IR-record / (IR-record + IR-record-calf))
      set IR-real food-available * split
      set IR-real-calf food-available * (1 - split)
    ]]


  ; reduce IR-record by IR
  set IR-record IR-record - IR-real
  if with-lact-calf? = true and wean-scale-fact != 1 [ set IR-record-calf IR-record-calf - IR-real-calf ]

  ; adjust intake rate based on fatness
  let max-SL (lgth * -0.39635) + 1.02347                                ; EQN 19: Calculated from Galatius & Kinze, unps.
  let over-mean-SL (storage-level - SL-mean) / (max-SL - SL-mean)       ; EQN 20: Check if storage levels exceed mean storage levels
  if over-mean-SL > 1 [ set over-mean-SL 1 ]
  let IR-SL-mod 1

  if over-mean-SL > 0 [                                                 ; If they do,
    let FC (-1 * satiation-c)
    set IR-SL-mod exp(over-mean-SL * FC)                                ; EQN 21: Calculate hunger modifier based on how much fatter they are than the average
    if IR-SL-mod < 0 [ set IR-SL-mod 0 ]                                ; Negative clamp
    set IR-real IR-real * IR-SL-mod                                     ; Adjust ingestion rate using the calculated modifier
  ]

  let food-eaten 0
  ifelse with-lact-calf? = FALSE or wean-scale-fact = 1 [set food-eaten IR-real][set food-eaten IR-real + IR-real-calf]

  ; Add to list to calculate daily ingestion
  set food-intake-list lput IR-real food-intake-list

  ; Remove eaten food from patches
  ask patch (item 0 (item 1 pos-list)) (item 1 (item 1 pos-list))  ; item 0 pos-list was the last added element, i.e. the current position
  [
      if food-level > 0 [
      set food-level food-level - food-eaten
      if ( food-level < 0.01 ) [ set food-level 0.01 ]
      ; if ( food-level > 0 and food-level <= 0.02 * maxU ) [ set pcolor 45 ]    ; turned off for efficiency
      ; if ( food-level > 0.02 and food-level <= 0.2 * maxU ) [ set pcolor 27 ]
      ; if ( food-level > 0.2 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 66 ]
      ; if ( food-level > 0.5 * maxU and food-level <= 0.9 * maxU ) [ set pcolor 63 ]
      ; if ( food-level > 0.9 * maxU and food-level < 1.1 * maxU ) [ set pcolor 61 ]
      ; if ( food-level > 1.1 * maxU ) [ set pcolor 61 ]
    ]
   ]

  set e-assim (IR-real * IR-to-EA * AE-food )                       ; EQN 22: Assimilate food of patch using the ingestion rate, the calibrated IR-real to EA conversion factor, and the assimilation efficiency of food
  set e-assim-calf (IR-real-calf * IR-to-EA * AE-food )             ; EQN 24: Calves assimilate food of patch using the ingestion rate, the calibrated IR-real to EA conversion factor, and the assimilation efficiency of food

  ; Debug energy intake
  if precision sim-day 3 = floor sim-day [if (debug = 11) [ if (who = 0) [ if food-available > 0 [
    print "Debugging energy intake: focal follow"
    print word "ID                 " who
    let preg 0
    ifelse m-preg > 0 [set preg "true"] [set preg "false"]
    print word "Pregnant?          " preg
    print word "With-lact-calf?    " with-lact-calf?
    print word "IR-struct-mass     " IR-struct-mass
    print (word "Repro supps       " " Preg: " preg-IR-sup " Lact: "Lact-IR-sup)
    print word "IR-timestep        " IR-timestep
    print word "Food-available     " food-available
    print word "IR-real            " IR-real
    print word "IR-real-calf       " IR-real-calf
    print word "over-mean-SL       " over-mean-SL
    print word "IR-SL-mod          " IR-SL-mod
    print word "food-eaten         " food-eaten
    print word "e-assim            " e-assim
    print word "e-assim-calf       " e-assim-calf
    print "END: Debugging energy intake: focal follow"
    print ""
  ]]]]

end


to calf-feed
  ; Calves over 3 mo also ingest food in patches to meet their metabolic needs
  let wean-scale-fact-calf (1.00 - wean-scale-fact) ; calf's portion of calf costs to cover

  ; Calf feeds and covers costs
  ifelse e-assim-calf >= (e-calf * wean-scale-fact-calf)                                                                                 ; Check if enough energy is available to cover total costs
  [ if (debug = 11) [if (who = 0) [ print word who "'s calf sufficiently fed - calf-feed" ]]
    set e-assim-calf e-assim-calf - (e-calf * wean-scale-fact-calf)                                                                      ; If so, deplete available energy to cover calf's share of costs
    set mass-struct-calf mass-struct-calf + (max-grow-calf * wean-scale-fact-calf)                                                       ; Calf grows maximally
    set v-blub-calf v-blub-calf + (((m-blub-calf * wean-scale-fact-calf)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))               ; Call adds storage energy to blubber
    if e-assim-calf > 0
    [ ifelse (v-blub-calf + ((e-assim-calf * DE-lip ) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub )) < (mass-calf * 0.436 / dens-blub) ; If extra energy available check if blubber volume is not at max (Mean + 1 SD percent blubber for calves in McLellan et al. 2002)
        [
        set v-blub-calf v-blub-calf + ((e-assim-calf * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                                  ; If not, add remaining energy to blubber volume
        set e-assim-calf 0
        ]
        [
        set e-assim-calf e-assim-calf - ((((mass-calf * 0.436 / dens-blub) - v-blub-calf)/ DE-lip ) * dens-blub * ED-lip)                ; If yes, set to max blubber volume
        set v-blub-calf (mass-calf * 0.436 / dens-blub)
        ]
    ]]

  [ ifelse e-assim-calf >= vital-costs-calf                                                                       ; If assimilated energy isn't enough to cover all of calf's portion of lactation costs, check if vital costs can be covered

    [ if (debug = 11) [if (who = 0) [  print word who "'s calf fed enough to cover calf portion of vital costs - calf-feed" ]]
      set e-assim-calf e-assim-calf - vital-costs-calf                                                            ; If storage isn't over mean levels then just cover vital lactation costs using e-assim-calf and split the remaining e-assim-calf to growth and blubber
      ifelse (e-assim-calf / 2) >= (m-growth-calf * wean-scale-fact-calf)                                         ; If half of remaining e-assim of the calf is more than enough to satisfy growth costs
      [ set e-assim-grow-calf (m-growth-calf * wean-scale-fact-calf) ]                                            ; Then set e-assim-calf for growth to the costs of growth and the remainder goes to calf storage
      [ set e-assim-grow-calf e-assim-calf / 2 ]                                                                  ; If not then e-assim-calf is split equally between growth and storage
            let e-assim-blub-calf e-assim-calf - e-assim-grow-calf                                                ; e-assim-calf for blubber is set to the remainder of e-assim-calf
            set e-assim-calf 0                                                                                    ; Set e-assim-calf to zero
            let growth-rate-calf (e-assim-grow-calf / ED-lean-mass) * DE-lean-mass                                ; Use available energy for growth
            set mass-struct-calf mass-struct-calf + growth-rate-calf
            set blub-calf ((e-assim-blub-calf / (ED-lip * perc-lip-blub) * DE-lip) / dens-blub)                   ; Calculate the blubber volume covered by available energy for blubber
            set v-blub-calf v-blub-calf + blub-calf                                                               ; Update calf blubber volume

    ]
    [ ifelse v-blub-calf >= (mass-calf * 0.313 / dens-blub)                                                                                            ; If e-assim-calf < vital-costs-calf check if calf is over mean percent blubber - 1 SD for calves in McLellan et al. 2002
          [ if (debug = 11) [if (who = 0) [ print word who "'s calf used blubber to cover costs - calf-feed" ]]
            let total-costs-calf e-calf - m-lact                                                                                                       ; Create a temp variable for the total calf costs
            let diff-calf-costs total-costs-calf - e-assim-calf                                                                                        ; Calculate the difference between this and the e-assim of the calf
            set v-blub-calf v-blub-calf - (((diff-calf-costs - (m-blub-calf * wean-scale-fact-calf))* DE-lip * perc-lip-blub) / (ED-lip * dens-blub )) ; Deplete blubber by difference
            set e-assim-calf 0                                                                                                                         ; Set e-assim-calf to zero
            set mass-struct-calf mass-struct-calf + (max-grow-calf * wean-scale-fact-calf)                                                             ; Calf grows maximally
           ]

          [                                                                                                     ; If e-assim-calf < vital-costs-calf, blubber volume is below mean, and energy is still needed pull energy for vital costs only from blubber
            if (debug = 11) [if (who = 0) [ print word who "'s calf used blubber to cover only vital costs - calf-feed" ]]
            let diff-calf-costs vital-costs-calf - e-assim-calf                                                 ; Find difference between vital costs and e-assim calf
            set v-blub-calf v-blub-calf - ((diff-calf-costs * DE-lip * perc-lip-blub) / (ED-lip * dens-blub))   ; Mobilize blubber to cover difference
            set e-assim-calf 0                                                                                  ; Set e-assim-calf to zero
           ]
      ]]

    check-calf-death
end

; ----------------- MAINTENANCE ----------------- ;
to maintenance
  ;;;;; COST CALCULATION ;;;;;
  set m-BMR (B0 * (weight ^ 0.75) * 1800)                                                       ; EQN 26: Calculate BMR using the current weight and B0 value, convert from watts to BMR timestep-1 by multiplying by 1800 seconds

  ;;;;; ALLOCATION ;;;;;
  ifelse e-assim >= m-BMR                                                                       ; Check if the energy assimilated is sufficient to cover BMR
    [ set e-assim e-assim - m-BMR ]                                                             ; If so reduce the available energy by the BMR cost
    [ set e-storage e-storage + e-assim                                                         ; If e-assim is not sufficient to cover BMR then add e-assim to storage and
      ifelse e-storage > m-BMR                                                                  ; Check if updated storage can cover BMR costs
    [ set e-storage e-storage - m-BMR                                                           ; If so, mobilize stored energy to do so
      set v-blub v-blub - (((m-BMR - e-assim)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub )) ; Reduce the amount of blubber volume by the BMR loss
    ]
    [                                                                                           ; If not, die
      set list-of-dead-age lput (floor age) list-of-dead-age
      set list-of-dead-day lput (floor sim-day) list-of-dead-day

      if (debug = 10) [ print word who " died of low body condition-Maintenance" ]
      die
    ]
      set e-assim 0                                                                             ; Set e-assim to zero
    ]


end


to calc-swim-speed ; calculates swimming speed for use in the thermoregulation and locomotion processes
  ifelse ( disp-type = 0 )
  [
    let swim-speed-pres-mov ((((10 ^ pres-logmov) * 100)) / 1800)
    let swim-speed-pres-log-trans ln(swim-speed-pres-mov * 100)
    set swim-speed swim-speed-pres-log-trans * 0.346                                    ; from movement data analysis-converting from swim-speed estimated on 30min to 6s scale from dead reckoned data
    if swim-speed <= 0.01 [set swim-speed 0.01]                                         ; negative clamp
  ]
  [
    set swim-speed (mean-disp-dist * 1000 * (1 /(30 * 60)))                             ; convert dispersal distance in km 30min-1 to m s-1
  ]
end

; ----- reporters to use in thermoregulation calculations ----- ;

; calculate heat transfer coefficient lower limit using a reporter - EQN 30
to-report calc-hC-low-lim[mean-blub-depth]
let hCll kB / mean-blub-depth                                                                                                          ; Hind and Gurney 1997
report hCll
end

; calculate forced convection scaling coefficient for each cone using a reporter - EQN 33
 to-report calc-for-convec-scal-coef[lengthcone]
let focsc (0.036) * (conductivity-w / ((lengthcone ^ (1 / 5)) * ( kin-visc-w ^ (4 / 5)))) * (Pr-w ^ (1 / 3)) * (swim-speed ^ (4 / 5))  ; Hind and Gurney 1997
report focsc
end

; ----- END: reporters to use in thermoregulation calculations ----- ;

; ----------------- THERMOREGULATION ----------------- ;
to thermoregulation
 ;;;;; COST CALCULATION ;;;;;
 let min-heat-gen-rate m-BMR + m-loco-ineff-pts                                                         ; EQN 28: Minimum heat generation rate incurred by maintenance and locomotion costs
 set m-loco-ineff-pts 0                                                                                 ; Clear locomotion inefficiency

 ;;; calculate for each cone ;;
 ; heat transfer coef lower limit
 (foreach (list (item 1 c2) (item 1 c3)) [                                                              ; Apply the heat transfer coefficient lower limit calculation to each of the cone mean blubber depths
    x -> let p calc-hC-low-lim(x)
    set hC-low-lim-coef-list lput p hC-low-lim-coef-list                                                ; Collect outputs in a list
    ])

 let heat-transfer-lower-lim mean hC-low-lim-coef-list                                                  ; Heat transfer lower limit as the average found for each cone
 set hC-low-lim-coef-list []                                                                            ; Clear list

 ; forced convection scaling coefficient
 (foreach (list (item 0 c2) (item 0 c3))  [                                                             ; Apply the forced convection scaling coefficient calculation to each of the cone lengths
   x -> let p calc-for-convec-scal-coef(x)
   set for-convec-scal-coef-list lput p for-convec-scal-coef-list                                       ; Collect outputs in a list
 ])

 let for-convec-scal-coef mean for-convec-scal-coef-list                                                ; Take the total forced convection scaling coefficient as the mean of the cone outputs
 set for-convec-scal-coef-list []                                                                       ; Clear list
 ;;; end of cone specific calculations ;;;

 let T-mb T-core - ((T-core - ([temp-w] of patch-here)) * 0.25)                                         ; EQN 32: Temperature of the blubber muscle interface relative to core and water temperature

 ; Calculating the skin temperature as a result of forced convection on the body - EQN 31
 let T-skin-1 ((11.4 * kB * lgth * T-mb) + (surface-area * temp-w * for-convec-scal-coef * ln(weight / mass-struct)))
 let T-skin-2 ((11.4 * kB * lgth) + (surface-area * for-convec-scal-coef * ln(weight / mass-struct)))
 let T-skin T-skin-1 / T-skin-2

; to calculate heat loss:
 let free-convec-scal-coef 0
 ifelse coef-thermal-expansion-w >= 0                                                                  ; Free convection scaling coefficient from Hind and Gurney 1997
  [set free-convec-scal-coef (0.555) * (conductivity-w / (surface-area / (2 * lgth)) ^ (1 / 4)) * (( gravity * coef-thermal-expansion-w * Pr-w ) / (kin-visc-w ^ 2))^(1 / 4)]
  [set free-convec-scal-coef (0.555) * (conductivity-w / (surface-area / (2 * lgth)) ^ (1 / 4)) * (( gravity * 0 * Pr-w ) / (kin-visc-w ^ 2))^(1 / 4)]
 let free-convec-heat-loss free-convec-scal-coef * (T-skin - temp-w) ^(3 / 2)                          ; Calculate the rate of heat loss from a free convection
 let for-convec-heat-loss (for-convec-scal-coef * (T-skin - temp-w))                                   ; Calculate the rate of heat loss from forced convection
 let rate-of-heat-loss free-convec-heat-loss + for-convec-heat-loss                                    ; Total rate of heat loss as the summation of heat loss from free and forced convection in W m-2

 let heat-trans-eff (min-heat-gen-rate / 1800) / (surface-area * (T-core - T-skin))                     ; EQN 34: Heat transfer efficiency
 let heat-gen-rate-real (surface-area * heat-transfer-lower-lim * (T-core - T-skin)) * 1800             ; EQN 29: Calculating Qcr or the realized heat generation rate
 let e-therm heat-gen-rate-real - min-heat-gen-rate                                                     ; Potential metabolic cost of thermo as the difference
 ifelse heat-trans-eff >= heat-transfer-lower-lim                                                       ; EQN 35: Check if the effective heat transfer coefficient is higher than the lower limit
 [ set m-thermo 0 ]                                                                                     ; If so, set thermoregulation costs to zero
 [ ifelse e-therm >= 0                                                                                  ; If not double check that the estimated metabolic cost of thermoregulation is greater than zero (negative check)
   [ set m-thermo e-therm ]                                                                             ; If it is, set thermo costs to the estimated metabolic costs
    [ set m-thermo 0 ] ]                                                                                ; If not, set to zero

 ;;;;; ALLOCATION ;;;;;
  ifelse e-assim >= m-thermo                                                                         ; check if the energy assimilated is sufficient to cover thermoregulatory costs
    [ set e-assim e-assim - m-thermo ]                                                               ; if so reduce the available energy by the thermo cost
    [ set e-storage e-storage + e-assim                                                              ; if e-assim is not sufficient to cover thermoregulation then add e-assim to storage and
      ifelse e-storage > m-thermo                                                                    ; check if updated storage can cover thermo costs
    [ set e-storage e-storage - m-thermo                                                             ; if so, mobilize stored energy to do so
      set v-blub v-blub - (((m-thermo - e-assim) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))  ; reduce the amount of blubber volume by the thermo loss
    ]
    [                                                                                                ; if not, die
      set list-of-dead-age lput (floor age) list-of-dead-age
      set list-of-dead-day lput (floor sim-day) list-of-dead-day

      if (debug = 10) [
        print word who " died of low body condition-Frozen"
      ]
      die
    ]
      set e-assim 0                                                                                          ; set e-assim to zero
    ]
end

; ----------------- LOCOMOTION ----------------- ;
to locomotion
  ;;;;; COST CALCULATION ;;;;;
  let aero-eff 0.13478 + 0.441 * ((swim-speed / 4.2) ^ 3) - 0.422 * ((swim-speed / 4.2) ^ 6)                    ; EQN 36: Aerobic efficiency calculation from Hind & Gurney 1997 modified for a 25 % max eff and 4.2 m s-1 max speed

  ; Calculate the Reynold's number
  let Re (lgth * swim-speed) / kin-visc-w                                                                       ; EQN 37: As in Fish 1998
  set Re log Re 10                                                                                              ; Convert to log scale

  ; Calculate CD using the log transformed relationship found for white sided dolphins in Tanaka et al 2019 - EQN 38
  let drag-coef (-0.113188 * Re + -1.535837)
  set drag-coef (10 ^ drag-coef)

  ; Calculate costs of locomotion
  set m-loco (lambda * density-w * surface-area * drag-coef * (swim-speed ^ 3)) / (2 * aero-eff * prop-eff)     ; EQN 39: As in Hind & Gurney 1997
  set m-loco m-loco * 1800                                                                                      ; Convert from watts to cost per timestep-1 by multiplying by 1800 seconds
  set m-loco-ineff-pts (1 - aero-eff) * m-loco                                                                  ; EQN 27: Save the waste heat generated this step to use in the thermo calculations next step

  ;;;;; ALLOCATION ;;;;;
  ifelse e-assim >= m-loco                                                                           ; Check if the energy assimilated is sufficient to cover locomotive costs
    [ set e-assim e-assim - m-loco ]                                                                 ; If so reduce the available energy by the locomotive cost
    [ set e-storage e-storage + e-assim                                                              ; If e-assim is not sufficient to cover locomotion then add e-assim to storage and
      ifelse e-storage >= m-loco                                                                     ; Check if updated storage can cover locomotive costs
    [ set e-storage e-storage - m-loco                                                               ; If so, mobilize stored energy to do so
      set v-blub v-blub - (((m-loco - e-assim)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))     ; Reduce the amount of blubber volume by the locomotive loss
    ]
    [ set e-storage (v-blub-min * dens-blub * ED-lip) ]
      set e-assim 0                                                                                  ; Set e-assim to zero
      ; For now this does not slow down swimming speed, but could in the future
    ]
end

; ----------------- REPRODUCTION ----------------- ;
to reproduction
   if with-lact-calf? = TRUE [
      lactation ]                                                                    ; Lactating females calculate lactation costs

  if pregnancy-status = 1 [
      pregnancy                                                                      ; Calculate pregnancy costs for pregnant females
      if round ds-mating =  t-gest [ give-birth ]]                                   ; Give birth
end

to porp-upd-pregnancy-status
  ; 0 (unable to mate, young/low energy); 1 (unable to mate, pregnant); 2 (ready to mate)

  ; Become ready to mate:
  if (pregnancy-status = 0 and age >= age-of-maturity ) [set pregnancy-status 2 ]         ; If of age set pregnancy status as ready to mate

  ; Mate:
  if (pregnancy-status = 2 and round (sim-day - 360 * (year - 1)) = mating-day) [         ; If ready to mate and mating season:
    let numT count(porps with [age >= age-of-maturity])                                   ; Total number of mature porpoises
    let numP round (numT * pregnancy-rate)                                                ; Use a temp variable to determine number of porpoises that should get pregnant
    let numR count porps with [v-blub >= v-blub-repro and age >= age-of-maturity]         ; Find total number above threshold that can get pregnant
    if preg-chance = 0 [set preg-chance random-float 1]

   ifelse numP <= numR                                                                    ; If there are more porpoises with a blubber volume higher than the reproductive threshold than the pregnancy rate
    [
      if (debug = 9 ) [ print "numR > numP" ]
      if (v-blub >= v-blub-repro) and (preg-chance <= pregnancy-rate) [
      set pregnancy-status 1                                                              ; These animals get pregnant
      set mass-f 0.000001                                                                 ; Set initial fetal mass "0"
      if (debug = 9 ) [ print word who " pregnant" ]
      set ds-mating 0 ]
    ]
    [
      ifelse (v-blub >= v-blub-repro)[                                                    ; If not enough animals are over reproductive threshold then have all animals over threshold
      set pregnancy-status 1                                                              ; Get pregnant
      set mass-f 0.000001                                                                 ; Set initial fetal mass "0"
      if (debug = 9 ) [ print word who " pregnant" ]   ; debug
      set ds-mating 0 ]
      [
      let upd-percent-pregnant ((numP - numR) / (numT - numR))                            ; Then some animals under threshold will get pregnant to reach 67% pregnancy rate
      if (preg-chance <= upd-percent-pregnant) [
      set pregnancy-status 1                                                              ; These animals get pregnant
      set mass-f 0.000001                                                                 ; Set initial fetal mass "0"
      if (debug = 9 ) [ print word who " pregnant" ]   ; debug
      set ds-mating 0 ]
      ]
  ]
]

  if pregnancy-status = 1 [ set ds-mating ds-mating + 1 ]                                 ; Update date in pregnancy period
  if with-lact-calf? [ set dsg-birth dsg-birth + 1 ]                                      ; Update date in the lactation period
end

to pregnancy
  ; Check if abortion/pregnancy occurred (without this sometimes a rare error will be produced)
  if mass-f <= 0 [
    set m-preg 0
    stop
  ]

  ;;;;; COST CALCULATION ;;;;;
  let max-grow-f (3 * (f-growth-c) ^ 3 ) * (( mass-f ^ (1 / 3) / f-growth-c) ^ 2)/ 48                                    ; EQN 41: Maximum fetal growth calculation ; / 48 for converting to time step
  set m-growth-f ((max-grow-f * percent-lip-f * ED-lip) / DE-lip) + ((max-grow-f * percent-pro-f * ED-pro) / DE-pro)     ; EQN 42: Calculate fetal tissue investment costs
  set e-heat-gest ((4400 * (max-mass-f ^ 1.2)) * 4184) / 14400                                                           ; EQN 43: Calculate heat of gestation in J per timestep
  set m-preg e-heat-gest +  m-growth-f                                                                                   ; EQN 40: Costs of pregnancy as the cost heat of gestation and fetal tissue investment

  ;;;;; ALLOCATION ;;;;;
  ifelse e-assim >= m-preg                                                                                               ; Check if enough energy is available to cover pregnancy costs
    [ set mass-f mass-f + max-grow-f                                                                                     ; If yes, then fetus grows maximally
      set e-assim e-assim - m-preg ]                                                                                     ; Deplete available energy by pregnancy costs
    [ set e-storage e-storage + e-assim                                                                                  ; If e-assim is not sufficient to cover pregnancy costs then add e-assim to storage and
      ifelse e-storage >= m-preg + e-repo-min                                                                            ; Check if updated storage can cover pregnancy costs and if it is more than necessary to continue reproducing
       [ set e-storage e-storage - m-preg                                                                                ; If so, mobilize stored energy to do so
         set mass-f mass-f + max-grow-f                                                                                  ; Fetus grows maximally
         set v-blub v-blub - (((m-preg - e-assim) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                     ; Reduce the amount of blubber volume by the pregnancy loss
         set e-assim 0 ]                                                                                                 ; Set e-assim to zero
    [ if ( debug = 9 ) [print word who " aborted calf"]                                                                  ; If not, then abort calf
         set pregnancy-status 0                                                                                          ; Set pregnancy status to 0 as too low energy to conceive
         set m-preg 0                                                                                                    ; Reset pregnancy costs to 0
         set abortion-count abortion-count + 1                                                                           ; Update abortion counter
         set n-calf-lost n-calf-lost + 1
         set mass-f 0 ] ]                                                                                                ; Reset fetal mass to default (0)
end

; birth
to give-birth
    set mass-calf mass-f                                               ; Set calf mass to the final fetal mass
    set mass-f 0                                                       ; Reset fetal mass
    set v-blub-calf (mass-calf * 0.375) / dens-blub                    ; Set the blubber volume of the calf to 37.5 % of mass (as in McLellan et al. 2002)
    set mass-struct-calf mass-calf * (1 - 0.375)                       ; Take the structural mass as the difference
    set pregnancy-status 2                                             ; It is ready to mate even though it has a very young calf
    set m-preg 0                                                       ; Reset pregnancy variables
    set e-heat-gest 0
    set m-growth-f 0
    set IR-record-calf 0
    set with-lact-calf? true
    let n random 2                                                     ; pick a random number to assign calf sex
    ifelse n = 1 [set sex-calf "female"] [set sex-calf "male"]
    set ds-mating -99                                                  ; reset days since mating
    set dsg-birth 0                                                    ; set days since giving birth as 0
    ask porps with [preg-chance > 0] [set preg-chance 0]
    if (debug = 9 ) [ print word who " with lact calf" ]
end

; lactation
to lactation
  ; Check if calf death occurred
  if mass-calf <= 0 [
    check-calf-death
    stop
    ]

  ;;;;; COST CALCULATION ;;;;;
  ;;;; calf maintenance - EQN 44 ;;;;
  set m-BMR-calf (B0 * (mass-calf ^ 0.75)) * 1800

  ; length calculations - EQN 45
  ifelse sex-calf = "female"
  [ let lgth-inf-f 116.3
    let lgth-0-f 89.2
    let lgth-k-f 3.3
    set lgth-calf ((lgth-inf-f * exp (ln(lgth-0-f / lgth-inf-f)* exp ((-1 * lgth-k-f)* (dsg-birth / 360))))/ 100)]                                                  ; Calculate length based on age for females; based on relationship from IDW dataset
  [ let lgth-inf-m 112.9
    let lgth-0-m 68.1
    let lgth-k-m 6.2
    set lgth-calf ((lgth-inf-m * exp (ln(lgth-0-m / lgth-inf-m)* exp ((-1 * lgth-k-m)* (dsg-birth / 360))))/ 100)]                                                  ; Calculate length based on age for males; based on relationship from IDW dataset

  ;;;; calf thermoregulatory costs - EQN 49 ;;;;
  let min-heat-gen-rate-calf m-BMR-calf                                                                                                                             ; Minimum heat generation rate incurred by maintenance costs; Assuming locomotion costs are negligable due to echelon swimming
  let for-convec-scal-coef-calf (0.036) * (conductivity-w / ((lgth-calf ^ (1 / 5)) *( kin-visc-w ^ (4 / 5)))) * (Pr-w ^ (1 / 3)) * (swim-speed ^ (4 / 5))           ; EQN 46: Scaling coefficient of forced convection from Hind and Gurney
  let storage-level-calf ((mass-calf - mass-struct-calf) / mass-calf)                                                                                               ; EQN 47: Calf storage level
  let dB-calf 6.21539 * storage-level-calf + 0.25118                                                                                                                ; EQN 48: Average blubber depth to storage level estimated for porpoise calves in the IDW dataset
  let heat-transfer-lower-lim-calf kB / dB-calf                                                                                                                     ; Calf heat transfer coefficient lower limit
  let surface-area-calf 0.093 * mass-calf ^ 0.57                                                                                                                    ; Calf surface area
  let T-mb-calf T-core - ((T-core - ([temp-w] of patch-here)) * 0.25)                                                                                               ; Temperature of the blubber muscle interface relative to core and water temperature
  let T-skin-calf-1 ((11.4 * kB * lgth-calf * T-mb-calf) + (surface-area-calf * ([temp-w] of patch-here) * for-convec-scal-coef-calf * ln(mass-calf / mass-struct-calf)))
  let T-skin-calf-2 ((11.4 * kB * lgth-calf) + (surface-area-calf * for-convec-scal-coef-calf * ln(mass-calf / mass-struct-calf)))
  let T-skin-calf T-skin-calf-1 / T-skin-calf-2                                                                                                                     ; Calculating the skin temperature as a result of forced convection on the body from Ryg et al. 1993 eqn 2

  let heat-trans-eff-calf (min-heat-gen-rate-calf / 1800) / (surface-area-calf * (T-core - T-skin-calf))                                                            ; Heat transfer efficiency
  let heat-gen-rate-real-calf (surface-area-calf * heat-transfer-lower-lim-calf * (T-core - T-skin-calf)) * 1800                                                    ; Calculating Qcr or the realized heat generation rate
  let e-therm-calf heat-gen-rate-real-calf - min-heat-gen-rate-calf                                                                                                 ; Potential metabolic cost of thermo as the difference
  ifelse heat-trans-eff-calf >= heat-transfer-lower-lim-calf                                                                                                        ; Check if the effective heat transfer coefficient is higher than the lower limit
  [ set m-thermo-calf 0 ]                                                                                                                                           ; If so, set thermoregulation costs to zero
  [ ifelse e-therm-calf >= 0                                                                                                                                        ; If not, double check that the estimated metabolic cost of thermoregulation is greater than zero (positive clamp)
    [ set m-thermo-calf e-therm-calf ]                                                                                                                              ; If it is, set thermo costs to the estimated metabolic costs
    [ set m-thermo-calf 0 ] ]                                                                                                                                       ; If not, set thermo to 0

  ;;;; calf growth costs - EQN 50 ;;;;
  let m-str-inf-c random-normal 15.43 0.59                                                                                                                          ; Von Bertalanffy fit from Galatius & Kinze, unps. dataset for calves < 0.9
  let m-str-k-c random-normal 20.95 4.07
  ifelse mass-struct-calf < m-str-inf-c
  [ set max-grow-calf (m-str-k-c / 17280) * ((m-str-inf-c ^ (1 / 3) * mass-struct-calf ^ (2 / 3)) - mass-struct-calf) ] [set max-grow-calf 0 ]                      ; Calculate max growth for calves
    set m-growth-calf (max-grow-calf * (ED-lean-mass + ED-lean-mass * (1 - DE-lean-mass)))                                                                          ; Calculate energy required for growth

  ;;;; calf blubber requirement costs - EQN 51 ;;;;
  ifelse v-blub-calf <  v-blub-calf-idl                                                                                                                             ; Check that calves blubber stores are sufficient
  [ set m-blub-calf (((v-blub-calf-idl - v-blub-calf) * dens-blub * perc-lip-blub * ED-lip) / DE-lip)]                                                              ; If blubber stores are not at ideal levels, calculate costs of bringing stores to that level
  [ set m-blub-calf 0 ]                                                                                                                                             ; If at sufficient levels, cost of blubber maintenance = 0

  ;;;; total calf costs - EQN 52 ;;;;
  set e-calf (m-BMR-calf + m-thermo-calf + m-growth-calf + m-blub-calf) / lact-eff                                                                                  ; calculate total calf costs
  if e-calf >= 337080.2 [                                                                                                                                           ; check if costs are over what mom can produce via milk in a timestep (as calculated from max values in Oftedal 1997)
    let m ((m-BMR-calf + m-thermo-calf + m-growth-calf) / lact-eff)                                                                                                 ; temp variable containing other costs
    set m-blub-calf (337080.2 - m) * lact-eff                                                                                                                       ; reduce m-blub-calf by other costs
    if m-blub-calf < 0 [set m-blub-calf 0]
  ]
  set e-calf (m-BMR-calf + m-thermo-calf + m-growth-calf + m-blub-calf) / lact-eff                                                                                  ; Recalculate total calf costs, should be less than or equal to 337080.2

  ; weaning scale factor
  if round dsg-birth <= (floor (0.375 * t-nurs))                                                                                                                    ; If calf is less than 3 months old, it is considered totally dependant ; MAYBE CHANGE THIS TO BE MASS DEPENDANT
  [ set m-lact e-calf ]                                                                                                                                             ; Calculate lactation costs for total coverage

  if (round dsg-birth > (floor (0.375 * t-nurs))) and (round dsg-birth <= t-nurs)                                                                                   ; If calf is less between 3 and 8 months old, it becomes less dependant with age
  [ let adj-date dsg-birth - (floor (0.375 * t-nurs))                                                                                                               ; Normalize scale
    set wean-scale-fact 1.00 - (adj-date / (0.625 * t-nurs))                                                                                                        ; Calculate the weaning scale factor for calves
    if wean-scale-fact < 0 [set wean-scale-fact 0 ]                                                                                                                 ; If for some reason the calf were to go over 8 mo set to zero
    set m-lact e-calf * wean-scale-fact ]                                                                                                                           ; Calculate lactation costs while considering weaning

  if round dsg-birth = t-nurs [ wean-calf ]                                                                                                                         ; If older than 8 months, fully wean calf

  ifelse wean-scale-fact < 1                                                                                                                                        ; Check if calf is being weaned
  [
  set vital-costs ((m-BMR-calf * wean-scale-fact) + (m-thermo-calf * wean-scale-fact))/ lact-eff                                                                    ; If so, set vital costs to maintenance and thermo costs offset by weaning scale factor
  set vital-costs-calf ((m-BMR-calf * (1.00 - wean-scale-fact)) + (m-thermo-calf * (1.00 - wean-scale-fact)))                                                       ; Set the calf's portion of vital costs to the remainder
  ]
  [
  set vital-costs (m-BMR-calf + m-thermo-calf)/ lact-eff                                                                                                            ; If not, set vital costs to maintenance and thermo costs
  ]

 if vital-costs > 337080.2 [                                                                                                                                        ; Check if vital costs are over what mom can produce via milk in a timestep (as calculated from maximum values in Oftedal 1997)
    set list-of-dead-age-calves lput (floor (dsg-birth / 360)) list-of-dead-age-calves
    set list-of-dead-day-calves lput (floor sim-day) list-of-dead-day-calves
    if (debug = 9) [print word who "'s calf died of high vital costs"]
    ; reset mother lactation variables
    set m-BMR-calf 0
    set lgth-calf 0
    set mass-calf 0
    set blub-calf 0
    set mass-struct-calf 0
    set m-thermo-calf 0
    set max-grow-calf 0
    set m-growth-calf 0
    set m-blub-calf 0
    set e-calf 0
    set m-lact 0
    set m-lact-real 0
    set vital-costs 0
    set vital-costs-calf 0
    set wean-scale-fact 1
    set with-lact-calf? false
    set dsg-birth -99
    set IR-record-calf 0
    set n-calf-lost n-calf-lost + 1
    stop
  ]

  ;;;;; ALLOCATION ;;;;;
  ifelse e-assim >= m-lact                                                                                                              ; Check if enough energy is available to cover lactation costs

  [ set e-assim e-assim - m-lact                                                                                                        ; If so, deplete available energy to cover lactation costs
    set mass-struct-calf mass-struct-calf + (max-grow-calf * wean-scale-fact)                                                           ; Calf grows maximally
    set v-blub-calf v-blub-calf + (((m-blub-calf * wean-scale-fact)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))
    set m-lact-real m-lact                                                                                                              ; All costs covered so realized cost = total lact cost
    ]

  [ ifelse e-assim >= vital-costs                                                                                                       ; If not enough to cover all lactation costs, check if vital costs can be covered

    [ ifelse v-blub >= v-blub-mean                                                                                                      ; If yes, first check if blubber volume is over mean values/animal is in good body condition
          [ set e-storage e-storage + e-assim                                                                                           ; If yes, then add e-assim to storage
            set e-storage e-storage - m-lact                                                                                            ; And then cover all costs using updated storage
            set v-blub v-blub - (((m-lact - e-assim) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                                 ; Reduce the amount of blubber volume by the lactation loss
            set mass-struct-calf mass-struct-calf + (max-grow-calf * wean-scale-fact)
            set v-blub-calf v-blub-calf + (((m-blub-calf * wean-scale-fact)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))
            set m-lact-real m-lact                                                                                                      ; All costs covered so realized cost = total lact cost
            set e-assim 0                                                                                                               ; Set e-assim to zero
      ]

          [ set m-lact-real e-assim                                                                                                     ; Only using e-assim to cover costs so realized cost = e-assim
            set e-assim e-assim - vital-costs                                                                                           ; If storage isn't over mean levels then just cover vital lactation costs using e-assim and split the remaining e-assim to growth and blubber
            ifelse (e-assim / 2) >= (m-growth-calf * wean-scale-fact)                                                                   ; If half of remaining e-assim is more than enough to satisfy growth costs
            [ set e-assim-grow-calf (m-growth-calf * wean-scale-fact) ]                                                                 ; Then set e-assim for growth to the costs of growth and the remainder goes to calf storage
            [ set e-assim-grow-calf e-assim / 2 ]                                                                                       ; If not then e-assim is split equally between growth and storage
            let e-assim-blub-calf e-assim - e-assim-grow-calf                                                                           ; e-assim for blubber is set to the remainder of e-assim
            set e-assim 0                                                                                                               ; Set e-assim to zero
            let growth-rate-calf (e-assim-grow-calf / ED-lean-mass) * DE-lean-mass                                                      ; Else, reduce growth rate and grow suboptimally; use only available energy for growth
            set mass-struct-calf mass-struct-calf +  growth-rate-calf
            set blub-calf ((e-assim-blub-calf / (ED-lip * perc-lip-blub)) * DE-lip)/ dens-blub                                          ; Calculate the blubber volume covered by available energy for blubber
            if e-assim-blub-calf >= (m-blub-calf * wean-scale-fact)                                                                     ; Check if added blubber would put the calf over ideal value
            [
              let blub-calf-adj (((e-assim-blub-calf - (m-blub-calf * wean-scale-fact))/ (ED-lip * perc-lip-blub)) * DE-lip)/ dens-blub ; If so, calculate how much over
              let blub-calf-dif blub-calf - blub-calf-adj                                                                               ; Calculate difference between initial blub-calf and adjusted blub-calf
              set e-assim (blub-calf-dif * dens-blub * perc-lip-blub * ED-lip)                                                          ; And add that energy back to e-assim
              set blub-calf blub-calf-adj                                                                                               ; Set blub-calf to adjusted value
            ]
            set v-blub-calf v-blub-calf + blub-calf ]]                                                                                  ; Update calf blubber volume


      [ ifelse v-blub >= v-blub-mean                                                                                                    ; If yes, first check if blubber volume is over mean values/animal is in good body condition
        [ set e-storage e-storage + e-assim                                                                                             ; If yes, then add e-assim to storage
          set e-storage e-storage - m-lact                                                                                              ; And then cover all costs using updated storage
          set v-blub v-blub - (((m-lact - e-assim) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                                   ; Reduce the amount of blubber volume by the lactation loss
          set mass-struct-calf mass-struct-calf + (max-grow-calf * wean-scale-fact)
          set v-blub-calf v-blub-calf + (((m-blub-calf * wean-scale-fact)* DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))
          set m-lact-real m-lact                                                                                                        ; All costs covered so realized cost = total lact cost
          set e-assim 0                                                                                                                 ; Set e-assim to zero
        ]

     [ ifelse v-blub >= v-blub-repro                                                                                                    ; If yes, first check if blubber volume is over mean values/animal is over reproductive threshold
        [ set e-storage e-storage + e-assim                                                                                             ; If e-assim is not sufficient to cover non-growth related lactation costs then add e-assim to storage and
          set e-storage e-storage - vital-costs                                                                                         ; And then cover vital costs using updated storage
          set v-blub v-blub - (((vital-costs - e-assim) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                              ; Reduce the amount of blubber volume by the vital lactation loss
          set m-lact-real vital-costs                                                                                                   ; Only vital costs covered so realized cost = vital costs
          set e-assim 0                                                                                                                 ; Set e-assim to zero
        ]
        [ set m-lact-real 0                                                                                                             ; No costs covered so realized cost = 0
          set v-blub-calf v-blub-calf - ((vital-costs * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))                                ; Have the calf use blubber to cover vital costs
          ]
  ]]]

  if wean-scale-fact < 1 [calf-feed]                                                                                                    ; If calf is being weaned trigger calf feeding procedure
end

to check-calf-death
  ; Calf should die if storage level falls below 0.05 - also has a probability of dying at higher storage levels in porp-upd-mortality
      let storage-level-calf 0
      if mass-calf > 0 [set storage-level-calf (v-blub-calf * dens-blub) / mass-calf ]
      if (storage-level-calf < 0.05) [
          set list-of-dead-age-calves lput (floor (dsg-birth / 360)) list-of-dead-age-calves
          set list-of-dead-day-calves lput (floor sim-day) list-of-dead-day-calves
          if (debug = 9) [print word who "'s calf died of low body condition-LACTATION"]
          ; reset mother lactation variables
          set m-BMR-calf 0
          set lgth-calf 0
          set mass-calf 0
          set blub-calf 0
          set mass-struct-calf 0
          set m-thermo-calf 0
          set max-grow-calf 0
          set m-growth-calf 0
          set m-blub-calf 0
          set e-calf 0
          set m-lact 0
          set m-lact-real 0
          set vital-costs 0
          set vital-costs-calf 0
          set wean-scale-fact 1
          set with-lact-calf? false
          set dsg-birth -99
          set IR-record-calf 0
          set n-calf-lost n-calf-lost + 1 ; for scenarios
        ]
end

to wean-calf
  let n-offspr 0                                                                                  ; Create temporary variable for the number of offspring
  if (sex-calf = "female" ) [ set n-offspr 1 ]                                                    ; Only female calves 'hatch'
  if (debug = 9 ) [
    let tmp word who " hatching "
    print word tmp n-offspr
  ]

  let age-calf dsg-birth / 360
  set dsg-birth -99

  hatch-porps  n-offspr [                                                   ; Create an independant porpoise agent representing weaned calf
    setxy random-xcor random-ycor
    set age age-calf                                                        ; Set porpoise age to calculated age
    set mass-struct mass-calf - (v-blub-calf * dens-blub)                   ; Set structural mass to difference between mass of calf and blubber mass
    set v-blub v-blub-calf                                                  ; Set blubber volume to calf blubber volume
    set weight mass-calf                                                    ; Set weight as calf mass
    set B0 11.13 + random-float 0.04 - random-float 0.04                    ; Normalization constant
    set lambda random-normal 0.25 0.016                                     ; Ratio of active to passive drag
    set DE-lip 0.74 + random-float 0.16                                     ; Deposition efficiency of lipid
    set DE-pro 0.43 + random-float 0.13                                     ; Deposition efficiency of protein
    set lgth-inf random-normal 158.12 4.68                                  ; Length asymptotic value
    set lgth-0 random-normal 94.82 1.69                                     ; Length initial value
    set lgth-k random-normal 0.41 0.06                                      ; Length k value
    set lgth lgth-calf                                                      ; Set length to calf length
    set surface-area 0.093 * weight ^ 0.57                                  ; Surface area
    set v-blub-min (weight * 0.05) / dens-blub                              ; set minimum blubber as 10% of weight
    set v-blub-repro (weight * repro-min-SL / dens-blub)                    ; set reproductive blubber volume threshold using minimum reproductive blubber volume
    set e-repo-min (v-blub-repro - v-blub-min) * dens-blub * ED-lip         ; calculate energy required for reproductive threshold
    set SL-mean ((lgth * -0.3059) + 0.7066)* IR-temp-mod                    ; mean storage level percentages from McLellan et al. 2002 for mature, calf, and immature porpoises
    set v-blub-mean (weight * SL-mean / dens-blub)
    set e-storage (v-blub - v-blub-min) * dens-blub * ED-lip                ; storage energy as the amount of energy stored in blubber over minimum threshold

  ; initialize reproductive costs
    set ds-mating -99                                                       ; Initial days since mating at -99
    set dsg-birth -99                                                       ; Initial days since giving birth at -99
    set pregnancy-status 0                                                  ; Unable to get pregnant, too young
    set with-lact-calf? false
    set m-preg 0
    set m-BMR-calf 0
    set lgth-calf 0
    set mass-calf 0
    set mass-struct-calf 0
    set m-thermo-calf 0
    set max-grow-calf 0
    set m-growth-calf 0
    set m-blub-calf 0
    set v-blub-calf 0
    set e-calf 0
    set m-lact 0
    set m-lact-real 0
    set vital-costs 0
    set vital-costs-calf 0
    set wean-scale-fact 1
    set IR-record-calf 0

  ; run upd-blubber-depths to initialize blubber values
    upd-blubber-depths
 ]

  ; reset mother lactation variables
     set m-BMR-calf 0
     set lgth-calf 0
     set mass-calf 0
     set mass-struct-calf 0
     set m-thermo-calf 0
     set m-growth-calf 0
     set m-blub-calf 0
     set e-calf 0
     set m-lact 0
     set m-lact-real 0
     set vital-costs 0
     set vital-costs-calf 0
     set wean-scale-fact 1
     set with-lact-calf? false
     set IR-record-calf 0
end

; ----------------- GROWTH ----------------- ;

to growth

   if mass-struct < m-str-inf                                                                                  ; If smaller than max female mass
      [
        ;;;;; COST CALCULATION ;;;;;
        set e-assim-grow e-assim / 2                                                                           ; Split available energy between growth and blubber

        set max-grow  (m-str-k / 17280) * ((m-str-inf ^ (1 / 3) * mass-struct ^ (2 / 3)) - mass-struct)        ; EQN 53: Maximum growth in a timestep ; 17280 adjuster to convert from annual to 30min basis

        ; calculate costs of growth
        set m-growth (max-grow * (ED-lean-mass + ED-lean-mass * (1 - DE-lean-mass)))                           ; EQN 54: Convert max growth of mass in kg to energy in joules per 30min
        if m-growth < 0 [set m-growth 0]                                                                       ; Negative clamp

        ;;;;; ALLOCATION ;;;;;
        ; EQN 57 within allocation possibilities
        ifelse e-assim-grow >= m-growth                                                                        ; Check if enough energy is available to cover completely
        [
          set e-assim e-assim - m-growth                                                                       ; If so, deplete available energy by the energy needed for growth
          set e-assim-grow 0                                                                                   ; Set e-assim for growth to zero
          set mass-struct mass-struct + max-grow                                                               ; Add growth mass to structural mass
          let growth-rate max-grow                                                                             ; Set growth rate to max growth rate
        ]
        [
          ifelse v-blub >= v-blub-mean                                                                         ; Check if blubber is in good condition
          [                                                                                                    ; If body condition > threshold grow maximally and pull extra needed energy from storage
            set e-storage e-storage + (e-assim-grow - m-growth)                                                ; Reduce e-storage by difference between e-assim for growth and energy needed for growth
            set e-assim e-assim -  e-assim-grow                                                                ; Reduce e-assim by e-assim for growth
            set v-blub v-blub - (((m-growth - e-assim-grow) * DE-lip * perc-lip-blub) / (ED-lip * dens-blub )) ; Reduce the amount of blubber volume by the growth difference
            set e-assim-grow 0                                                                                 ; Set e-assim for growth to zero
            set mass-struct mass-struct + max-grow                                                             ; Grow maximally
            let growth-rate max-grow                                                                           ; Set growth rate to max growth rate
           ]
           [
              let growth-rate (e-assim-grow / ED-lean-mass) * DE-lean-mass                                     ; Else, reduce growth rate and grow suboptimally; use only available energy for growth
              set mass-struct mass-struct + growth-rate                                                        ; Grow by adjusted growth rate
              set e-assim e-assim -  e-assim-grow                                                              ; Reduce e-assim by e-assim for growth
              set e-assim-grow 0                                                                               ; Set e-assim for growth to zero
           ]
        ]
      ]

  let lgth-t ((lgth-inf * exp (ln(lgth-0 / lgth-inf)* exp((- lgth-k)* age)))/ 100)                             ; EQN 58: Update length based on age
  if lgth-t > lgth [set lgth lgth-t]

 end


; ----------------- STORAGE ----------------- ;
to upd-blubber-depths
  ; Called once per day to update the site specific blubber depths of porps based on their storage levels and lengths - EQN 60
  ; See TRACE section 2.7.2.3.3 & 3.2.2 for details

  let SL-over-L storage-level / (lgth * 100)           ; Storage level adjusted by length

  let average-dB v-blub / (surface-area * 10000)       ; EQN 61: Blubber depth if it were to be distributed evenly over body surface

  ;;; Site name:                                             Short abb:    Long abb:
  ; axillary dorsal site ----------------------------------- AxD --------- axillary-D
  ; axillary lateral site ---------------------------------- AxL --------- axillary-L
  ; axillary ventral site ---------------------------------- AxV --------- axillary-V
  ; cranial insertion of the dorsal fin dorsal site -------- CrD --------- CrIDF-D
  ; cranial insertion of the lateral fin dorsal site ------- CrL --------- CrIDF-L
  ; cranial insertion of the ventral fin dorsal site ------- CrV --------- CrIDF-V
  ; caudal insertion of the dorsal fin dorsal site --------- CaD --------- CaIDF-D
  ; caudal insertion of the lateral fin dorsal site -------- CaL --------- CaIDF-L
  ; caudal insertion of the ventral fin dorsal site -------- CaV --------- CaIDF-V

  ; Site specific blubber depth slope and intercepts
  let q-bl-AxD (random-normal 31.35 13.46)
  let q-bl-AxL (random-normal -23.86 11.53)
  let q-bl-AxV (random-normal -45.16 14.84)
  let q-bl-CrD (random-normal 47.36 13.80)
  let q-bl-CrL (random-normal 9.80 10.21)
  let q-bl-CrV (random-normal -29.33 14.27)
  let q-bl-CaD (random-normal -14.62 16.12)
  let q-bl-CaL (random-normal 16.39 12.87)
  let q-bl-CaV (random-normal 13.56 14.36)
  let b-bl-AxD (random-normal 1.42 0.05)
  let b-bl-AxL (random-normal 1.58 0.04)
  let b-bl-AxV (random-normal 1.76 0.05)
  let b-bl-CrD (random-normal 1.53 0.05)
  let b-bl-CrL (random-normal 1.51 0.04)
  let b-bl-CrV (random-normal 1.66 0.05)
  let b-bl-CaD (random-normal 1.68 0.06)
  let b-bl-CaL (random-normal 1.46 0.05)
  let b-bl-CaV (random-normal 1.42 0.05)

  ; Calculate scaling factor value based on storage level
  let Ax-D-a q-bl-AxD * SL-over-L + b-bl-AxD
  let Ax-L-a q-bl-AxL * SL-over-L + b-bl-AxL
  let Ax-V-a q-bl-AxV * SL-over-L + b-bl-AxV
  let CrIDF-D-a q-bl-CrD * SL-over-L + b-bl-CrD
  let CrIDF-L-a q-bl-CrL * SL-over-L + b-bl-CrL
  let CrIDF-V-a q-bl-CrV * SL-over-L + b-bl-CrV
  let CaIDF-D-a q-bl-CaD * SL-over-L + b-bl-CaD
  let CaIDF-L-a q-bl-CaL * SL-over-L + b-bl-CaL
  let CaIDF-V-a q-bl-CaV * SL-over-L + b-bl-CaV

  ; Blubber depth calculation
  set axillary-D Ax-D-a * average-dB
  set axillary-L Ax-L-a * average-dB
  set axillary-V Ax-V-a * average-dB
  set CrIDF-D CrIDF-D-a * average-dB
  set CrIDF-L CrIDF-L-a * average-dB
  set CrIDF-V CrIDF-V-a * average-dB
  set CaIDF-D CaIDF-D-a * average-dB
  set CaIDF-L CaIDF-L-a * average-dB
  set CaIDF-V CaIDF-V-a * average-dB

  ; Cone outputs as a list
  set c2 (list (lgth * .20) (((axillary-D + axillary-L + axillary-V + CrIDF-D + CrIDF-L + CrIDF-L)/ 6)/ 100)) ; Cone 2
  set c3 (list (lgth * .21) (((CrIDF-D + CrIDF-L + CrIDF-L + CaIDF-D + CaIDF-L + CaIDF-L)/ 6)/ 100))          ; Cone 3

 end

;    --------------------
;    |      UPDATE      |
;    --------------------

to storage

  ; if any excess energy remains after allocation, save it to storage
  if e-assim >= 0 [                                                                        ; if any remaining energy
    let add-Blub ((e-assim * DE-lip * perc-lip-blub) / (ED-lip * dens-blub ))              ; EQN 59: convert assimilated energy to blubber volume
    set v-blub v-blub + add-Blub                                                           ; update blubber volume
    set e-assim 0]                                                                         ; set e-assim to 0


  ; update blubber depths once per day
  if ((remainder time-step 48) = 0) [
    set storage-level-sum 0                                                                ; reset daily ; CHANGED BY CARA
    upd-blubber-depths
  ]

end

to upd-state-variables

  set m-tot m-BMR + m-loco + m-thermo + m-growth + m-preg + m-lact-real        ; EQN 25: Calculate total expended energy

  set weight mass-struct + (v-blub * dens-blub)                                ; EQN 65: Update weight

  set storage-level (weight - mass-struct) / weight                            ; EQN 66: Storage level update
  set storage-level-sum storage-level-sum + storage-level

  ifelse pregnancy-status != 1                                                 ; EQN 67: Calculate surface area - Worthy and Edwards 1990
  [set surface-area 0.093 * weight ^ 0.57]
  [set surface-area 0.093 * (weight + (2 * mass-f)) ^ 0.57]

  set SL-mean ((lgth * -0.3059) + 0.7066) * IR-temp-mod                        ; EQN 68: Mean storage level percentages from McLellan et al. 2002 for mature, calf, and immature porpoises
  set v-blub-mean (weight * SL-mean / dens-blub)                               ; EQN 69: Mean blubber volume
  set v-blub-min (weight * 0.05) / dens-blub                                   ; EQN 70: Update minimum blubber volume
  set v-blub-repro (weight * repro-min-SL / dens-blub)                         ; EQN 71: Update reproductive blubber volume threshold, mean - 2*SD from McLellan et al. 2002

  set e-storage (v-blub - v-blub-min) * dens-blub * ED-lip                     ; EQN 72: Update storage energy using the updated blubber volume and minimum blubber volume
  set e-repo-min (v-blub-repro - v-blub-min) * dens-blub * ED-lip              ; EQN 73: Update reproductive energy threshold

  if with-lact-calf? = TRUE
  [
    set mass-calf mass-struct-calf + (v-blub-calf * dens-blub)                 ; EQN 74: Update calf mass
    set v-blub-calf-idl (mass-calf * calf-idl-SL / dens-blub)                  ; EQN 75: Mean percent blubber for calves in McLellan et al. 2002
  ]

end

;    -----------------------
;    |      SCENARIOS      |
;    -----------------------


to scenario-outputs
  set daily-pop-list lput (count porps) daily-pop-list
  let current-ships ss-ships with [start-ts <= time-step and end-ts >= time-step]
  if any? current-ships
  [ ask current-ships
    [ set fourty-km-list lput (count porps in-radius 100) fourty-km-list
      if any? porps in-radius 100 [ set storage-level-deter-list lput (mean [storage-level] of porps in-radius 100) storage-level-deter-list ]
      if any? porps with [with-lact-calf? = TRUE] in-radius 100 [set mean-calf-mass-deter-list lput ( mean [mass-calf] of porps with [with-lact-calf? = TRUE]) mean-calf-mass-deter-list]
      if any? porps with [age < age-of-maturity] in-radius 100 [ set mean-juv-mass-deter-list lput ( mean [weight] of porps with [age < age-of-maturity]) mean-juv-mass-deter-list]
       ]
   ]

  if sim-day = 5040 [ set n-calf-lost 0 ]

  if sim-day = 5400 [
    set tot-calf-loss-out n-calf-lost
    set fourty-km-out fourty-km-list
    set storage-level-out storage-level-deter-list
    set mean-calf-mass-out mean-calf-mass-deter-list
    set mean-juv-mass-out mean-juv-mass-deter-list
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
398
12
1006
1021
-1
-1
1.0
1
7
1
1
1
0
1
1
1
0
599
0
999
0
0
0
ticks
30.0

BUTTON
28
18
140
51
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
153
18
260
51
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
1014
287
1106
332
disp-var
disp-var
"bathymetry" "maxent-level" "food-prob" "food-level" "blocks" "temperature" "salinity"
3

MONITOR
28
64
91
109
NIL
time-step
0
1
11

MONITOR
101
64
164
109
NIL
sim-day
1
1
11

CHOOSER
1013
521
1205
566
debug
debug
"profile" 0 1 2 3 4 5 6 7 8 9 10 11 12
1

CHOOSER
1013
476
1205
521
model
model
0 1 2 3 4
4

SLIDER
1014
46
1206
79
n-porps
n-porps
1
500
200.0
1
1
NIL
HORIZONTAL

CHOOSER
1014
78
1206
123
area
area
"Kattegat" "Homogeneous"
0

BUTTON
1106
286
1206
332
Update disp-var
landsc-display
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
173
64
236
109
NIL
year
17
1
11

MONITOR
247
64
310
109
NIL
quarter
1
1
11

INPUTBOX
1014
123
1206
183
max-sim-day
14400.0
1
0
Number

PLOT
31
140
387
376
population
time (days)
NIL
0.0
900.0
0.0
3000.0
true
true
"" ""
PENS
"N x10" 1.0 0 -10899396 true "" ""
"E x100" 1.0 0 -14070903 true "" ""
"total food" 1.0 0 -2674135 true "" ""

MONITOR
323
319
382
364
N porps
count porps
17
1
11

PLOT
30
379
209
499
porpoise-storage-level
NIL
NIL
0.05
0.6
0.0
20.0
true
false
"" ""
PENS
"default" 0.025 1 -14070903 true "" ""

CHOOSER
1014
243
1206
288
write-data
write-data
"off" "daily" "monthly" "yearly" "one-porp" "one-porp-energy-debug" "all-porps-energy-debug"
6

INPUTBOX
1014
183
1206
243
output-name
outfile
1
0
String

TEXTBOX
1015
458
1197
476
Model selection and testing
12
0.0
1

CHOOSER
1014
331
1206
376
wind-farms
wind-farms
"off" "Pot_Krieger" "Pot_St_Middelgr" "Rodsand-I" "Rodsand-II" "Samsoe" "Sprogoe" "User-def" "All" "Line" "Potential"
0

SWITCH
1014
376
1206
409
incl-ships?
incl-ships?
0
1
-1000

SLIDER
1012
924
1204
957
deterrence-coeff
deterrence-coeff
0
10.0
0.7
0.1
1
NIL
HORIZONTAL

TEXTBOX
1014
838
1164
856
Noise avoidance param.
12
0.0
1

SLIDER
1012
956
1204
989
std-deterrence-dist
std-deterrence-dist
200
40000
40000.0
20
1
m
HORIZONTAL

SLIDER
1012
988
1204
1021
deter-time
deter-time
1
10
5.0
1
1
time steps
HORIZONTAL

MONITOR
323
270
382
315
N lact
count porps with [ with-lact-calf? = true ]
0
1
11

MONITOR
321
64
384
109
NIL
month
1
1
11

TEXTBOX
1016
769
1201
799
Mortality param.
12
0.0
1

TEXTBOX
330
201
389
219
in patches
11
0.0
1

SLIDER
1013
788
1205
821
bycatch-prob
bycatch-prob
0
0.2
0.0
0.001
1
y-1
HORIZONTAL

SWITCH
1014
409
1206
442
pile-driving?
pile-driving?
0
1
-1000

MONITOR
324
221
382
266
N preg
count porps with [pregnancy-status = 1]
0
1
11

PLOT
31
769
387
896
field-metabolic-rates
NIL
NIL
0.0
360.0
12.0
30.0
true
false
"" ""
PENS
"NonRepro" 1.0 0 -2674135 true "" ""
"Repro" 1.0 0 -8630108 true "" ""

PLOT
208
379
386
499
age-distribution
year class
NIL
0.0
25.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -15040220 true "" ""

PLOT
208
498
386
618
mass-length
NIL
NIL
1.0
1.7
0.0
4.0
true
false
"" ""
PENS
"default" 1.0 2 -14730904 true "" ""

PLOT
30
498
210
618
blubber-mass
repro class
NIL
0.0
1.0
0.0
25.0
true
true
"" ""
PENS
"Calves" 1.0 2 -13791810 true "" ""
"Juvs" 1.0 2 -13840069 true "" ""
"Nonlact" 1.0 2 -5825686 true "" ""
"Lact" 1.0 2 -8630108 true "" ""

TEXTBOX
32
117
182
136
Model outputs:
15
0.0
1

TEXTBOX
1015
746
1165
765
Parameters:
15
0.0
1

SLIDER
1012
892
1204
925
resp-level
resp-level
100
200
165.0
1
1
dBpp
HORIZONTAL

SLIDER
1012
860
1204
893
source-level
source-level
150
300
262.0
1
1
dBpp
HORIZONTAL

CHOOSER
1013
597
1205
642
scenario
scenario
0 1 2 3 4 5 6 7 8 9
0

CHOOSER
1013
642
1205
687
sc-month
sc-month
1 2 3 4 5 6 7 8 9 10 11 12
1

PLOT
31
642
387
770
blubber-depths
NIL
NIL
0.0
360.0
1.5
3.5
true
false
"" ""
PENS
"axillary-D" 1.0 0 -7500403 true "" ""
"axillary-L" 1.0 0 -2674135 true "" ""
"axillary-V" 1.0 0 -955883 true "" ""
"CrIDF-D" 1.0 0 -6459832 true "" ""
"CrIDF-L" 1.0 0 -1184463 true "" ""
"CrIDF-V" 1.0 0 -10899396 true "" ""
"CaIDF-D" 1.0 0 -13840069 true "" ""
"CaIDF-L" 1.0 0 -14835848 true "" ""
"CaIDF-V" 1.0 0 -11221820 true "" ""

PLOT
31
894
387
1021
energy-intake
NIL
NIL
0.0
360.0
7.5
40.0
true
false
"" ""
PENS
"Adults" 1.0 0 -2674135 true "" ""
"Repro" 1.0 0 -8630108 true "" ""

TEXTBOX
1015
580
1165
598
Scenario selectors
12
0.0
1

CHOOSER
1013
687
1205
732
sc-year
sc-year
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
19

TEXTBOX
1017
20
1167
38
Model run selections:
14
0.0
1

TEXTBOX
33
622
221
640
Seasonally varying outputs:
13
0.0
1

BUTTON
272
18
385
51
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## ODD FOR THE UPDATED HARBOUR PORPOISE POPULATION MODEL

Model and description updated by Cara Alyse Gallagher 2019-01-10

The ODD protocol presented here includes elements from the model descriptions of four different harbor porpoise models. The original model (M1) was developed to explicitly model the movements of porpoises in Danish waters (Nabe-Nielsen et al. 2013). The second model (M2) extended M1 by incorporating simplistic energetics and modelled the effects of anthropogenic noise and bycatch on the porpoise population in the inner Danish waters (Nabe-Nielsen et al. 2014). The third model (M3) was developed as a further extension of the model into the North Sea, with area-specific large-scale movements and more realistic behavioral effects of disturbances (Nabe-Nielsen et al. 2018); it is important to note that the noise and transmission loss modelling from M3 was incorporated into the model presented here, but the movements were retained from M1 and M2, as these correspond to animals in the inner Danish waters. M4 represents the current version of the model presented here. 


## PURPOSE

The purpose of the model is to simulate the population dynamics of harbor porpoises, Phocoena phocoena, in the Inner Danish Waters (IDW), an environment with large seasonal fluctuations in temperature and salinity, and to investigate how energy budgets can be used to evaluate the seasonal variation in impacts of noise in this species, in particular seismic surveys. This model is an extension to a previously published IBM of harbor porpoises in the IDW (Nabe-Nielsen et al. 2013, 2014) as it represents the individual energetics of harbor porpoises in a more realistic manner. Our approach is based on first principles and explicitly accounts for variations in energy requirements among individual animals. 

## STATE VARIABLES AND SCALES

The model is composed of three kinds of entities: porpoise agents (called “porps”), disturbance agents (here ships), and landscape grid cells. As in Nabe-Nielsen et al. (2014), porp agents are characterized by their location, speed, movement direction, age, age of maturity, storage level, pregnancy status and lactation status (ODD Table. 2.1.). Here we also include the state variables length, mass, blubber mass, and body condition. Porps are ‘super-individuals’ (Scheffer et al. 1995), hereafter referred to as individuals, that represent several female porpoises with characteristics aligning with those of porpoises in the IDW population. Ships are represented by their location, speed, movement direction, and noise level, as in Nabe-Nielsen et al. (2014). The simulated landscape is also conserved from Nabe-Nielsen et al. (2014) and covers a 240 km x 400 km area of the IDW centered around the islands Funen and Zealand (see ODD Fig. 2.1). The landscape is composed of 600 x 1000 grid cells, each covering 400m x 400m.

Each landscape grid cell is characterized by its location, average water depth, distance to land, food level and maximum food level (as in Nabe-Nielsen et al. 2014) (Table. 2.2.); in this study we introduce the state variables location-specific water temperature and salinity to allow for spatiotemporally specific calculation of the thermophysical parameters involved in metabolic calculations. Grid cells are divided into land (52.1%), water without food (47.1%), and randomly distributed food patches (4572 individual patches) (0.76%), as in Nabe-Nielsen et al. (2014). The maximum food level for each food patch was estimated by assuming that its food level correlated with the probability of porpoise presence as determined by maximum entropy (MaxEnt) estimates of porpoise distribution for animals tagged with satellite transmitters (Edrén et al. 2010, Nabe-Nielsen et al. 2014). There is no food outside of food patches. The porps’ large-scale movements are guided by the maximum season-specific food level in different 40km x 40km squares, called “blocks”. 

The model proceeds in discrete time steps of 30-minutes and typically runs for 20 years. Each month consists of 30 days and a year of 360 days. Metabolic calculations are in units of energy per unit time (J 30min-1).

## PROCESS OVERVIEW AND SCHEDULING

M2: At the beginning of each timestep porps move and respond to noise by adjusting their movements if disturbed (Movement) (see Fig. 2.2 for visual process overview). The movement model includes memory of favorable patches to which movement is oriented when moving in unfavorable locations. Thereby, dynamic homeranges emerge, reflecting the distribution of resources and conspecifics.  If an animal’s storage level has decreased for 3 days, it switches to large-scale movements. These movements aren’t affected by noise.

After moving, porps ingest and expend energy on their metabolic costs (Energy budget). Energy is ingested and assimilated from prey patches if available. Assimilated and stored energy are used to allocate to metabolic processes in order of importance to survival, following Sibly et al. (2013). Any energy remaining after fulfilling costs is stored as blubber, affecting the storage level of the animal. 
  
M1: The food levels in the patches are dynamically changing. When a patch is visited by a porp its food level drops corresponding to the amount consumed by the porp, but afterwards it increases logistically in daily steps until reaching the maximum food level for the patch. M2: Whereas the logistic food growth rate is kept constant, the maximum food level in a patch depends on the time of the year and on where it is located. The season-specific maximum food level was obtained from maximum entropy (MaxEnt) estimates of where porpoises are most likely to occur (Edrén et al., 2010). These were calculated from the number of observations of satellite-tracked animals in areas with different environmental conditions.

Once every day porps may die (with probabilities that depend on their energy level and age), mate, become pregnant, give birth, or wean a calf (depending on the time of the year and reproductive status; submodel Life-history). 

Environmental parameters, such as water temperature and salinity, are updated monthly (on the 1st day of each month) and correspond to location-specific monthly averages from the Copernicus database for 2007-2016 (http://marine.copernicus.eu). The thermophysical properties of seawater are then determined for each patch using lookup tables specific to each temperature and salinity combination. The formulae used to determine the values used in the lookup tables can be found in Table 3.3 in TRACE section 3. 
M2: If day of year is 60, 150, 210 or 300: Season-specific MaxEnt data are updated. Between days 300 and day 60 (i.e. during winter): the values are rescaled to lie in the range 0–1. For other seasons: values are rescaled to achieve the same potential food level for the entire landscape (to ensure that the sum of the maximum food levels for all patches does not change among seasons). If the day of year is 1: each porp calculates a new mating date.
The porp agents’ blubber mass depends on their movements, morphometrics (size and shape), and reproductive state. State variables related to morphometrics (V_bl, m, SA, SL) and 〖IR〗_max are updated at the end of every time step (30 minutes; Table 2.1.), while other variables are updated asynchronously; i.e. immediately after execution of their related process. Body site-specific blubber depth, d_(B,site), and age are updated once per day. 


## DESIGN CONCEPTS

Basic principles:
This model is built on the assumptions that the porpoise population in the IDW is food limited and that noise affects porpoises by displacing them from specific regions, hence reducing their food intake. Energy budgets are based on physiological principles (Sibly et al. 2012, Sibly et al. 2013) and mathematical models of energy expenditure (Hind and Gurney 1997, Gallagher et al. 2018). This model builds on the ability of animals to make behavioral and reproductive decisions based on memory (Van Moorter et al. 2009) and current body condition, such as when to move to a previously visited food patch or move to a new region or when to abort or abandon calves. 

Emergence:
Population size and demographics emerge from adaptive decisions made by porps related to their movements and energetic state, which have the potential to ultimately improve individual fitness. The realized growth rates of individuals, age at first successful reproduction, and body condition are emergent properties of the model and affect the size and demographics of the population. While the use of a pregnancy rate imposes a maximum cap on the annual reproductive effort, the actual number of calves born in a year emerges from the current environmental and physiological state of females. Age-related and seasonal fluctuations in body condition emerge from variations in foraging efficiency and environmental conditions experienced by porps. While the potential maximum energy taken in in a timestep is modified using a water-temperature based modifier, the total energy consumed and stored is based on the resources encountered and the energetic costs experienced by animals. The relationship between mortality probability and storage level is imposed but an animal’s probability of dying at any point in time is driven by the storage level, which emerges from the animal’s foraging success. The spatial distribution and home ranges of porps in the model emerge from the way individuals move in response to changing energy levels. 

Objectives:
Animals attempt to maximize their fitness by maximizing energy intake through memory-based movement and by allocating available energy to processes in order of necessity. If processes necessary to survival are not covered by energy intake, stored energy will be used to cover them. As storage levels decrease, the animal has an increasing probability of dying. Processes not critical to immediate survival can be left unmet with some trade-off occurring (Sibly et al. 2013), such as reduced growth or calf loss. 

Prediction:
Animals base their prediction of how much food they can gather in different food patches on their previous visits to those patches, and this influences their fine-scale movements. Their memory, though, is limited and, after some time, patches visited in the distant past are forgotten.  Animals are assumed to have knowledge about which parts of the landscape have the highest food levels at different times of the year. This permits them to move to one of the 40 km × 40 km blocks with the highest potential food level when using large-scale movements.

Sensing:
Porps can sense their storage levels and make decisions specific to their current state, e.g. move to higher quality areas, reduce reproductive investment, or use energy converted from stores. They can sense water depth and temperature. Sensing water depth allows them to turn away from coast lines. Porps can also sense noise levels from ships and turbines, and this affects their movements.

Interaction:
The porps only interact indirectly by foraging on the same food patches. Once they begin to take solid foods (beginning at 3mo) dependent calves foraging in the same patches also compete with their mothers. 

Stochasticity:
Fine-scale movements, individual morphometrics, energetic inputs, mortality, and dates related to reproduction (e.g. mating date, calving date, and date of weaning) are associated with stochasticity in the model (see Table 3.1). In all these cases, the variation represented is considered important, as it may have consequences, but its mechanistic basis is unknown or considered not relevant for the purpose of the model. The probability of surviving increases with increasing energy levels.

Collectives:
Mothers and dependent calves in the model act as collectives. They move and forage as one unit. Super-individuals represent approximately 100 individual porpoises, assessed using the latest population counts (SCANS III, 2016), but their energetics are modelled to resemble the energetics of a single animal; they were introduced for run-time reasons and to represent densities which are typical of the IDW. 

Observation:
The number of animals, their storage levels, and the total amount of food in the landscape are recorded daily. Age-class distributions and age-specific mortalities are recorded yearly. 


## INITIALISATION

The model is by default initialized by creating 200 randomly distributed female porp individuals, representing 20,000 real porpoises. Their age-class distribution corresponds to that of stranded and bycaught animals (Lockyer & Kinze 2003) and 67% of the adults are pregnant. No porps are initialized with lactating calves. Their initial morphometric state variables (length, mass, storage level, and site-specific blubber depth) correspond to those of by-caught female porpoises from the Inner Danish and nearby waters (unpublished data from C. Kinze & A. Galatius (IDW porpoise dataset)) (TRACE Section 3 Table 3.3). Mating, calving, and weaning periods of the individuals were set to occur at day 225 (±20), day 182 (±20), and day 60 (±20), respectively (mean ±SD; random normal variables; Lockyer & Kinze 2003, Lockyer 2003, Nabe-Nielsen et al. 2014). Simulations were initiated on 1 January and food levels in food patches were initialized as the location specific maximum food levels for that date. Location-specific environmental parameters were taken as the average found for the month of January from the Copernicus dataset for the years 2007–2016. If running simulations, survey ships are initialized and the XY positions of the “buoys” used as turning points for the vessels to follow the survey route are loaded at setup. Parameters related to noise and disturbance can be chosen on the user interface, with defaults found in the Movement model parameter table (Table 2.3.).

The procedures involved in initialization are: 
•	Setup
•	Setup_Global_Parameters
•	Landsc_Setup
•	Porps_Setup
•	Update_Block_Values
•	Ships_Setup

In the Setup procedure, which is called by the UI “Setup” button, the world is cleared and completely reset, aside from the interface slider values. In the Setup_Global_Parameters procedure all global parameter values, e.g. constants such as the energy densities of tissues, gestation and lactation period length, etc., are set. In Landsc_Setup procedure, all environmental parameters are set up and maps are imported. In Porps_Setup, the model first sets up porp parameters and then assigns their state variables. After Porps_Setup, the Update_Block_Values procedure is run to initialize the mean MaxEnt values for blocks for each quarter. If a scenario is selected using the Netlogo “chooser” scenario, and a scenario month and scenario order (for Scenario 2) are selected using the “sc_month” and “sc_order” choosers, the survey ships are then setup using the Ships_Setup procedure. 


## INPUT DATA
Seven types of background maps are used in the model. The food probability (a binomial 0 or 1 value) is determined using the file ‘patches.asc’. The bathymetry and distance to coast for patches is set using the files ‘bathy.asc’ and ‘disttocoast.asc’, respectively, and block number is determined using the file ‘blocks.asc’. The maximum food level of a food patch was obtained from seasonal MaxEnt estimates, as in Nabe-Nielsen et al. (2014) (ascii files: ‘quarter’1-4’.asc’). Environmental conditions, like water temperature and salinity, were obtained from the Copernicus database Baltic Sea Reanalysis Product for years 2007-2016 and updated once a month (ascii files: ‘TempData_Month’X’.asc’ & ‘SalData_Month’X’.asc’) (see TRACE Data Evaluation for more information on temperature and salinity maps). All ascii maps used have a spatial resolution equal to the patch size of 400m × 400m, however, as the dataset used to create the temperature and salinity maps has a spatial resolution of 4km × 4km, patches falling within each 4km × 4km contain identical temperature and salinity values.

For each patch, the thermophysical properties of seawater (i.e. water density, thermal conductivity, coefficient of thermal expansion, specific heat capacity, dynamic viscosity, Prandtl number, and kinematic viscosity) are pulled from lookup tables (loaded at setup and stored as global matrices) once a month following the temperature and salinity map update using the location-specific water temperature and salinity values (in the ‘thermophysical properties of SW tables’ folder: ‘DensityTable.csv’, ‘ThermalConductivity.csv’, ‘SpecificHeat.csv’, ‘CoefOfThermalExpansionTable.csv’, & ‘DynamicVis.csv’). The Prandtl number, Pr, and kinematic viscosity of water, ν, are then determined using these values (see TRACE Section 3 Table 3.4 for details).

Simulations include details for simulating disturbance from seismic surveys. These details are provided in a txt file with the location identifier of the survey route, x,y coordinates of navigation buoys, source level, and start and end timestep of the survey. Noise is emitted for all timesteps where a survey vessel is present. The source level of all surveys was taken as 262 dB re 1 µPa peak-peak (a standard metric for noise levels in water; Thomspson et al. 2013, Kyhn et al. 2019). The threshold where a behavioral reaction occurred was set at 165 dB re 1 µPa peak-peak corresponding to the minimum distance recorded to result in a behavioral reaction in wild porpoises responding to a commercial seismic survey in Thompson et al. (2013). See the Movement model for details on how noise from seismic survey ship agents was represented in the model. 


## SUBMODELS

The model is composed of two distinct submodels, the 1) Movement model and the 2) Energy Budget model

1) The movement model used here was developed in Nabe-Nielsen et al. (2013) and employed in Nabe-Nielsen et al. (2014) and is based on memory-based adaptive foraging behavior. This model was informed and calibrated using movement track data from porpoises in the inner Danish waters (Sveegaard et al. 2011). The Movement model allows for porps to make adaptive movement decisions based on spatial memory and their current energetic state and to sense and respond to anthropogenic disturbances.

2) The energy budget submodel developed here serves as an application and extension of the equation-based model in Gallagher et al. (2018) within an individual-based populations model. It is based on established general principles of physiological ecology. The model combines elements of the methods presented in Sibly et al. (2013) and Hind and Gurney (1997) to mechanistically model the energetic requirements of marine homeotherms. Sibly et al. (2013) serves as the foundation for the model framework and the allocation process follows this approach while the methods employed to model hydrodynamic and thermal processes follow those of Hind & Gurney (1997). This model links the spatiotemporal variability in environmental conditions to the energy balance, survival probabilities, and reproductive success of cetaceans to be used as a tool in assessing the population level impacts of anthropogenic disturbance with high spatiotemporal resolution.   




## REFERENCES

Agricultural Research Council, 1980. The Nutrient Requirements of Ruminant Livestock. Slough, Berks.

Anderson, S.S. and Fedak, M.A., 1987. Grey seal, Halichoerus grypus, energetics: females invest more in male offspring. Journal of Zoology, 211(4), pp.667-679.

Barham, R.J., 2016. Omø South Nearshore A/S: Underwater noise. Orbicon A/S. Technical Report from Orbicon & Subacoustech Environmental Ltd No. OS-TR-003. URL: https://ens.dk/sites/ens.dk/files/Vindenergi/os-tr-003_underwater_noise.pdf

Beltran, R.S., Testa, J.W. and Burns, J.M., 2017. An agent-based bioenergetics model for predicting impacts of environmental change on a top marine predator, the Weddell seal. Ecological Modelling, 351, pp.36-50.

Blanchet, M., Nance, T., Ast, C., Wahlberg, M. and Acquarone, M., 2008. First case of a monitored pregnancy of a harbour porpoise (Phocoena phocoena) under human care. Aquatic Mammals, 34(1), p.9.

Blaxter, K., 1989. Energy metabolism in animals and man. CUP Archive.

Boult, V.L., Quaife, T., Fishlock, V., Moss, C.J., Lee, P.C. and Sibly, R.M., 2018. Individual-based modelling of elephant population dynamics using remote sensing to estimate food availability. Ecological modelling, 387, pp.187-195.

Brody, S. 1968. Bioenergetics and growth. Hafner Publishing Company Inc. New York, NY.

Brown, J.H., Gillooly, J.F., Allen, A.P., Savage, V.M. and West, G.B. 2004. Toward a metabolic theory of ecology. Ecology, 85, 1771–1789.

Camphuysen, C. J., and Krop, A. 2011. Maternal care, calf-training and site fidelity in a wild harbour porpoise in the North Sea. Lutra 54:123-126.

Edrén, S., Wisz, M.S., Teilmann, J., Dietz, R. and Söderkvist, J., 2010. Modelling spatial patterns in harbour porpoise satellite telemetry data using maximum entropy. Ecography, 33(4), pp.698-708.

Fish, F.E., 1993. Power output and propulsive efficiency of swimming bottlenose dolphins (Tursiops truncatus). Journal of Experimental Biology, 185(1), pp.179-193.

Fish, F.E. 1996. Transitions from drag-based to lift-based propulsion in mammalian swimming. American Zoologist, 36(6), pp.628-641.

Fish, F.E., 1998. Comparative kinematics and hydrodynamics of odontocete cetaceans: morphological and ecological correlates with swimming performance. Journal of Experimental Biology, 201(20), pp.2867-2877.

Gallagher, C.A., Stern, S.J. and Hines, E., 2018.The metabolic cost of swimming and reproduction in harbor porpoises (Phocoena phocoena) as predicted by a bioenergetic model. Marine Mammal Science. Doi: https://doi.org/10.1111/mms.12487

George, J., 2009. Growth, morphology and energetics of bowhead whales, Balaena mysticetus, University of Alaska Fairbanks, Fairbanks. PhD thesis.  URL: http://hdl.handle.net/11122/9031

Gibbs, C.L. and Gibson, W.R., 1972. Energy production of rat soleus muscle. American Journal of Physiology-Legacy Content, 223(4), pp.864-871.

Gol'din, P. E. 2004. Growth and body size of the harbour porpoise, Phocoena phocoena (Cetacea, Phocoenidae), in the Sea of Azov and the Black Sea.

Hind, A.T. and Gurney, W.S., 1997. The metabolic cost of swimming in marine homeotherms. Journal of Experimental Biology, 200(3), pp.531-542. 

Kastelein, R.A., Helder-Hoek, L. and Jennings, N., 2018. Seasonal changes in food consumption, respiration rate, and body condition of a male harbor porpoise (Phocoena phocoena). Aquatic Mammals, 44(1), pp.76-91.

Kleiber, M., 1975. The fire of life: an introduction to animal energetics. Rev. ed. Pub 

Krieger RE, Malabar, Florida.

Koopman, H.N., 1998. Topographical distribution of the blubber of harbor porpoises (Phocoena phocoena). Journal of Mammalogy, 79(1), pp.260-270.

Koopman, H.N., Pabst, D.A., Mclellan, W.A., Dillaman, R.M. and Read, A.J., 2002. Changes in blubber distribution and morphology associated with starvation in the harbor porpoise (Phocoena phocoena): evidence for regional differences in blubber structure and function. Physiological and Biochemical Zoology, 75(5), pp.498-512.

Kriete, B., 1994. Bioenergetics in the killer whale, orcinus orca, University of British Columbia, Vacncouver. PhD thesis.  Doi: 10.14288/1.0088104

Kyhn, L.A., Wisniewska, D.M., Beedholm, K., Tougaard, J., Simon, M., Mosbech, A. and Madsen, P.T., 2019. Basin-wide contributions to the underwater soundscape by multiple seismic surveys with implications for marine mammals in Baffin Bay, Greenland. Marine pollution bulletin, 138, pp.474-490.

Lockyer, C., 1991. Body composition of the sperm whale, Physeter catodon, with special reference to the possible functions of fat depots. Marine Research Institute.

Lockyer, C. 1995. Investigation of aspects of the life history of the harbour porpoise, Phocoena phocoena, in British waters. Report of the International Whaling Commission, pp.189-199.

Lockyer, C., 2003. Harbour porpoises (Phocoena phocoena) in the North Atlantic: Biological parameters. NAMMCO Scientific Publications, 5, pp.71-89.

Lockyer, C. and Kinze, C., 2003. Status, ecology and life history of harbour porpoise (Phocoena phocoena), in Danish waters. NAMMCO Scientific Publications, 5, pp.143-175.

Malavear, M.Y.G., 2002. Modeling the energetics of Steller sea lions (Eumetopias jubatus) along the Oregon coast, Oregon State University, Corvallis. MSc thesis. URL: https://ir.library.oregonstate.edu/concern/graduate_projects/wd3761431

McLellan, W.A., Koopman, H.N., Rommel, S.A., Read, A.J., Potter, C.W., Nicolas, J.R., Westgate, A.J. and Pabst, D.A., 2002. Ontogenetic allometry and body composition of harbour porpoises (Phocoena phocoena, L.) from the western North Atlantic. Journal of Zoology, 257(4), pp.457-471.

Nabe-Nielsen, J., Tougaard, J., Teilmann, J., Lucke, K., Forchhammer, M., 2013. How a simple adaptive foraging strategy can lead to emergent home ranges and increased food intake. Oikos 122, 1307–1316.

Nabe-Nielsen, J., Sibly, R.M., Tougaard, J., Teilmann, J. and Sveegaard, S., 2014. Effects of noise and by-catch on a Danish harbour porpoise population. Ecological Modelling, 272, pp.242-251. 

Nabe-Nielsen, J., van Beest, F.M., Grimm, V., Sibly, R.M., Teilmann, J. and Thompson, P.M., 2018. Predicting the impacts of anthropogenic disturbances on marine populations. Conservation Letters, 11(5), p.e12563. 

NIRAS, Rambøll, and DHI, 2015. Underwater noise and marine mammals. Technical Report from 
Energinet.dk No. 15/06201-1. Rev. nr. 4. URL: https://naturstyrelsen.dk/media/162585/underwater-noise-and-marine-mammals_2392023_rev4.pdf

Oftedal, O. T. 1997. Lactation in whales and dolphins: evidence of divergence between baleen and toothed species. Journal of Mammary Gland Biology and Neoplasia 2(3):205-230.

Otani, S., Naito, Y., Kato, A. and Kawamura, A., 2001. Oxygen consumption and swim speed of the harbor porpoise Phocoena phocoena. Fisheries science, 67(5), pp.894-898.

Parry, D. A. 1949. The structure of whale blubber, and a discussion of its thermal properties. Journal of Cell Science 3(9):13-25.

Pinheiro, J., Bates, D., DebRoy, S. and Sarkar, D., 2019. R Core Team (2017) nlme: linear and nonlinear mixed effects models. R package version 3.1-139. Computer software: Retrieved from https://CRAN. R-project. org/package= nlme

Prentice, A.M. and Prentice, A., 1988. Energy costs of lactation. Annual review of nutrition, 8(1), pp.63-79.

Pullar, J.D. and Webster, A.J.F., 1977. The energy cost of fat and protein deposition in the rat. British Journal of Nutrition, 37(3), pp.355-363.

R Core Team., 2019. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. http://www.R-project.org/

Read, A.J. and Hohn, A.A., 1995. Life in the fast lane: the life history of harbor porpoises from the Gulf of Maine. Marine Mammal Science, 11(4), pp.423-440.

Rechsteiner, E.U., Rosen, D.A. and Trites, A.W., 2013. Energy requirements of Pacific white-sided dolphins (Lagenorhynchus obliquidens) as predicted by a bioenergetic model. Journal of Mammalogy, 94(4), pp.820-832.

Rojano-Doñate, L., McDonald, B.I., Wisniewska, D.M., Johnson, M., Teilmann, J., Wahlberg, M., Højer-Kristensen, J. and Madsen, P.T., 2018. High field metabolic rates of wild harbour porpoises. Journal of experimental biology, 221(23), p.jeb185827.

Ryg, M., Smith, T.G. and Øritsland, N.A., 1988. Thermal significance of the topographical distribution of blubber in ringed seals (Phoca hispida). Canadian Journal of Fisheries and Aquatic Sciences, 45(6), pp.985-992.

Ryg, M., Lydersen, C., Knutsen, L.Ø., Bjørge, A., Smith, T.G. and Øritsland, N.A., 1993. Scaling of insulation in seals and whales. Journal of Zoology, 230(2), pp.193-206.

SCANS-III, 2016. Small Cetaceans in the European Atlantic and the North Sea (SCANS-III). Final Report. University of St Andrews, UK. http://biology.st-andrews.ac.uk/scans3/.

Scheffer, M., Baveco, J.M., DeAngelis, D.L., Rose, K.A. and van Nes, E., 1995. Super-individuals a simple solution for modelling large populations on an individual basis. Ecological modelling, 80(2-3), pp.161-170.

Sibly, R.M., Brown, J.H. and Kodric-Brown, A. eds., 2012. Metabolic ecology: a scaling approach. John Wiley & Sons.

Sibly, R.M., Grimm, V., Martin, B.T., Johnston, A.S., Kułakowska, K., Topping, C.J., 

Calow, P., Nabe‐Nielsen, J., Thorbek, P. and DeAngelis, D.L., 2013. Representing the acquisition and use of energy by individuals in agent‐based models of animal populations. Methods in Ecology and Evolution, 4(2), pp.151-161. 

Sørensen, T.B. and Kinze, C.C., 1994. Reproduction and reproductive seasonality in Danish harbour porpoises, Phocoena phocoena. Ophelia, 39(3), pp.159-176.

Sveegaard, S., Teilmann, J., Tougaard, J., Dietz, R., Mouritsen, K.N., Desportes, G. and 

Siebert, U., 2011. High‐density areas for harbor porpoises (Phocoena phocoena) identified by satellite tracking. Marine Mammal Science, 27(1), pp.230-246.

Tanaka, H., Li, G., Uchida, Y., Nakamura, M., Ikeda, T., Liu, H. 2019. Measurement of time-varying kinematics of a dolphin in burst accelerating swimming. PLoS ONE 14(1): e0210860. https://doi.org/10.1371/journal.pone.0210860

Thompson, P.M., Brookes, K.L., Graham, I.M., Barton, T.R., Needham, K., Bradbury, G. & Merchant, N.D. 2013. Short-term disturbance by a commercial two-dimensional seismic survey does not lead to long-term displacement of harbour porpoises. - Proceedings of the Royal Society B-Biological Sciences 280:8.

Turchin, P., 1998. Quantitative analysis of movement. Sinauer.

Urick, R.J., 1983. Propagation of sound in the sea: Transmission loss, I and II. Principles of underwater sound. McGraw-Hill Inc, New York, pp.177-179.

Van Moorter, B., Visscher, D., Benhamou, S., Börger, L., Boyce, M.S. and Gaillard, J.M., 2009. Memory keeps you at home: a mechanistic model for home range emergence. Oikos, 118(5), pp.641-652.

Worthy, G. A., and E. F. Edwards. 1990. Morphometric and biochemical factors affecting heat loss in a small temperate cetacean (Phocoena phocoena) and a small tropical cetacean (Stenella attenuata). Physiological Zoology 63(2):432-442.

Worthy, G.A., 1991. Insulation and thermal balance of fasting harp and grey seal pups. Comparative Biochemistry and Physiology Part A: Physiology, 100(4), pp.845-851.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

boat top
true
14
Polygon -7500403 true false 150 1 137 18 123 46 110 87 102 150 106 208 114 258 123 286 175 287 183 258 193 209 198 150 191 87 178 46 163 17
Rectangle -16777216 false true 129 92 170 178
Rectangle -16777216 false true 120 63 180 93
Rectangle -7500403 true false 133 89 165 165
Polygon -11221820 true false 150 60 105 105 150 90 195 105
Polygon -16777216 false true 150 60 105 105 150 90 195 105
Rectangle -16777216 false true 135 178 165 262
Polygon -16777216 false true 134 262 144 286 158 286 166 262
Line -16777216 true 129 149 171 149
Line -16777216 true 166 262 188 252
Line -16777216 true 134 262 112 252
Line -16777216 true 150 2 149 62

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="110210 Test disp-dist" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="52416"/>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-sim-110210&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;one-porp&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;bathymetry&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110317  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.5"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
      <value value="400"/>
      <value value="500"/>
      <value value="600"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110322a  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.5"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110322b  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="16"/>
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="200"/>
      <value value="300"/>
      <value value="400"/>
      <value value="500"/>
      <value value="600"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404a Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404b Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404c Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404d Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110405 Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="345600"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406a Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406b Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406c Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406d Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406e Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.25"/>
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406f Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.5"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110411 Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412a Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412b Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412c Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412d Calibrate e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412e Calibrate e-use (REFERENCE)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412f Calibrate e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference model" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + ships + windf" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + pot wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + pot wf constr" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference model + extra juv mort" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-025" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-050" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-100" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-000" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111007 - Calibrate adult mort m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="m-mort-prob-const" first="0.3" step="0.1" last="0.7"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-5" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-6" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-7" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-8" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-9" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference, no bycatch" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-10" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-11" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + pot wf" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + pot wf + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.01" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.017" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.05" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.5" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference, no bycatch, but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + pot wf + all ships but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + pot wf but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference, no bycatch, but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + pot wf + all ships but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + pot wf but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference, no bycatch, but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf + all ships but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference, no bycatch, but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf + all ships but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120419 New reference + all wind farms + bycatch=0.1" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120419 New reference + all wind farms + bycatch=0.2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120503 New reference, no bycatch, but rU=0.005" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms + pot wf + all ships but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms + all ships but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + pot wf + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Potential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + bycatch=0.05" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + all wind farms + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + potential wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Potential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=4" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=3" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=1" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="211218 ingestion calibration" repetitions="2" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="259200"/>
    <exitCondition>count porps &gt; 350 or count porps &lt; 100</exitCondition>
    <metric>cal_param</metric>
    <metric>count porps</metric>
    <metric>mean [storage-level] of porps</metric>
    <metric>mean [caidf_l] of porps</metric>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="460"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="5400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR_Coef">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;temperature&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="IR2EA" first="110000000" step="2000000" last="150000000"/>
    <enumeratedValueSet variable="pregnancy_rate">
      <value value="0.68"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_min_blub_perc">
      <value value="0.156"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fatness_coef">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110119 ingestion calibration" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="259200"/>
    <exitCondition>count porps &gt; 350 or count porps &lt; 100</exitCondition>
    <metric>cal_param</metric>
    <metric>count porps</metric>
    <metric>mean [storage-level] of porps</metric>
    <metric>mean [caidf_l] of porps</metric>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="460"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="5400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR_Coef">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;temperature&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="IR2EA" first="110000000" step="1000000" last="120000000"/>
    <enumeratedValueSet variable="pregnancy_rate">
      <value value="0.68"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_min_blub_perc">
      <value value="0.156"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fatness_coef">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120119 ingestion calibration" repetitions="2" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="259200"/>
    <exitCondition>count porps &gt; 350 or count porps &lt; 100</exitCondition>
    <metric>cal_param</metric>
    <metric>count porps</metric>
    <metric>mean [storage-level] of porps</metric>
    <metric>mean [caidf_l] of porps</metric>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="460"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="5400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR_Coef">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;temperature&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="IR2EA" first="105000000" step="1000000" last="112000000"/>
    <enumeratedValueSet variable="pregnancy_rate">
      <value value="0.68"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_min_blub_perc">
      <value value="0.156"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fatness_coef">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="010519 Scenario 2 Test" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="346560"/>
    <metric>thirty_km_list</metric>
    <metric>storage_level_list</metric>
    <metric>mean_calf_mass_list</metric>
    <metric>mean_juv_mass_list</metric>
    <metric>mean_calf_mass_deter_list</metric>
    <metric>mean_juv_mass_deter_list</metric>
    <metric>tot_calf_loss</metric>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc_month">
      <value value="3"/>
      <value value="6"/>
      <value value="9"/>
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;Test&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR2EA">
      <value value="185000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="21000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="5400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR_Coef">
      <value value="4.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satiation_c">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resp_level">
      <value value="176"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_min_blub_perc">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1.285"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="source_level">
      <value value="0"/>
      <value value="262"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="030519_Scen3_Jan" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="345600"/>
    <metric>daily_pop_list</metric>
    <metric>thirty_km</metric>
    <metric>storage_level_deter</metric>
    <metric>mean_calf_mass_deter</metric>
    <metric>mean_juv_mass_deter</metric>
    <metric>tot_calf_loss</metric>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc_month">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;Test&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR2EA">
      <value value="185000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="21000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="7200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="IR_Coef">
      <value value="4.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satiation_c">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resp_level">
      <value value="176"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_min_blub_perc">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1.285"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="source_level">
      <value value="0"/>
      <value value="262"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
