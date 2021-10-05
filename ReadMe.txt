Pendulum.osim: OpenSim Model file
Pendulum2D.scone: Scone file used for optimizations.
PendulumController.lua: Basic controller file, with optimized open-loop control input
PendulumController_NoCocon.lua: Controller file with open-loop control input equal to 0.01
PendulumController_CoconFix.lua: Controller file with open-loop input fixed to predefined value
PendulumController_internalnoise.lua: Controller file where noise is added to control signal instead of as a moment
PendulumController_internalnoise_nococon.lua: Same as before, open-loop input equal to 0.01
PendulumControllerSeparate: Controller file where separate parameters are optimized for muscle 1 and muscle 2
PendulumMeasureBarrier: Objective used throughout
