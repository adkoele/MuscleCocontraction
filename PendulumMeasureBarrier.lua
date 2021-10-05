-- script to calculate objective for the muscle-driven  inverted pendulum
-- Ton van den Bogert 4/21/2021

function init( model )
	
	-- cost function weights
	if not ( scone.Wtrack and scone.Weffort and scone.RMS) then
		error( "Must set Wtrack and Weffort parameters!" )
	end	
	Wtrack = tonumber(scone.Wtrack)
	Weffort = tonumber(scone.Weffort )
	RMS = tonumber(scone.RMS)
	max_dur = tonumber(scone.dur)
	
	-- Find the dof and muscles
	angle = model:find_dof( "q1")
	muscle1 = model:actuator( 1 )
	muscle2 = model:actuator( 2 )
	
	-- Preallocate a table in which we will record the simulation data
	local max_duration = 100
	local update_interval = 0.001;
	local size = math.tointeger(max_duration/update_interval) + 1
	sim = { angle={}, u1={}, u2={} }
	for i = 1,size do
		sim.angle[ i ] = 0.0
		sim.u1[ i ] = 0.0
		sim.u2[ i ] = 0.0
	end
	simdataptr = 0
	
end

function update(model)
	
	-- Record angle and the two muscle controls
	simdataptr = simdataptr + 1
	sim.angle[simdataptr]   = angle:position()
	sim.u1[simdataptr] = muscle1:input()
	sim.u2[simdataptr] = muscle2:input()
	
	
	--local s = string.format('PendumumMeasure update at t=%8.3f', 	model:time() )
	--scone.info(s)
	
end

function result(model)
	
	duration = model:time() 	
	
	sum1 = 0.0	
	sum2 = 0.0
	for i=1, simdataptr do
		sum1 = sum1 + (sim.angle[i]-3.14159266 )^2
		sum2 = sum2 + sim.u1[i]^2 + sim.u2[i]^2
	end
	tracking_cost = sum1/simdataptr
	effort_cost = sum2/simdataptr/2
	
	if tracking_cost > (RMS/180*3.14159266)^2 then
		cost = 100+Wtrack*tracking_cost + Weffort*effort_cost
	else
		cost = Weffort*effort_cost
	end
	
	local s = string.format('RMSangle(deg):%7.3f RMSu:%7.3f', 180/math.pi*math.sqrt(tracking_cost), math.sqrt(effort_cost) )
	scone.info(s)
	local s = string.format('Tracking term:%9.3f Effort term:%7.3f', Wtrack*tracking_cost, Weffort*effort_cost )
	scone.info(s)
	local s = string.format('simdataptr: %d  angle: %f', simdataptr, sim.angle[simdataptr] )
	scone.info(s)	
	local s = string.format('Duration: %f', duration )
	scone.info(s)
	
	-- add a penalty if duration < max_duration
	max_duration = max_dur
	if duration < max_duration then
		duration_penalty = 100*(max_duration - duration)
		cost = cost + duration_penalty
		local s = string.format('Duration penalty: %f', duration_penalty )
		scone.info(s)
	end
	
	return cost
end

