-- SCONE script to control the inverted pendulum
-- adapted from Tutorial 6b - Script - Balance Device

function init( model, par )
	
	-- find the model elements we need in the controller
	target_body = model:find_body( "rod1" )
	angle = model:find_dof( "q1")
	muscle1 = model:actuator( 1 )
	muscle2 = model:actuator( 2 )
	
	-- get controller parameters that are being optimized
	if ( scone.Kp and scone.Kd ) then
		Kp = par:create_from_string( "Kp", scone.Kp )
		Kd = par:create_from_string( "Kd", scone.Kd )
		delay = par:create_from_string( "delay", scone.delay )
	else
		error( "Must set u0, Kp, Kd, delay parameters!" )
	end	
	
	-- perturbation parameters
	if not ( scone.perturb_interval and scone.perturb_amplitude ) then
		error( "Must set perturb_interval and perturb_amplitude parameters!" )
	end	
	perturb_interval = tonumber(scone.perturb_interval)
	perturb_amplitude = tonumber(scone.perturb_amplitude )
	perturb_seed = tonumber(scone.perturb_seed)
	
	-- global variables to keep track of perturbation, and the last time that it was changed
	perturb_time = -perturb_interval  -- to force a perturbation at t=0
	perturb_moment = 0.0
	
	-- make sure we use the same perturbation signal in each evaluation
	math.randomseed( perturb_seed )
	
	-- Preallocate a buffer for delayed feedback
	buffer_size = 200
	buffer = { t = {}, ang={}, vel={} }
	for i = 1,buffer_size do
		buffer.t[ i ] = 0
		buffer.ang[ i ] = 0
		buffer.vel[ i ] = 0
	end
	buffer_ptr = 0
	
end

function update( model )
	
	local t = model:time()
	local ang = 0
	local vel = 0
	
	--store angle and velocity in buffer
	buffer_ptr = buffer_ptr + 1
	if buffer_ptr > buffer_size then
		buffer_ptr = 1
	end
	buffer.t[buffer_ptr]   = t
	buffer.ang[buffer_ptr] = angle:position()
	buffer.vel[buffer_ptr] = angle:velocity()
	
	--find the two samples in the buffer that bracket the time t-delay
	local i1=0
	local i2=0
	local d1min=1e5
	local d2min=1e5
	for i=1,buffer_size do
		d1 = t-delay - buffer.t[i]
		if d1>=0 and d1<d1min then
			d1min=d1
			i1=i
		end
		d2 = buffer.t[i] - (t-delay)
		if d2>=0 and d2<d2min then
			d2min=d2
			i2=i
		end
	end
	--scone.debug( "i1=" .. i1 .. "  i2=" .. i2 .."  d1min=" .. d1min .. "  d2min=" .. d2min)
	
	--interpolate between i1 and i2 to extract the delayed feedback
	if i1==0 then         -- this means that we have no feedback yet
		ang = 3.14159266 
		vel = 0
	elseif i1==i2 then
		ang = buffer.ang[i1]
		vel = buffer.vel[i1]
	else
		ang = ( d1min*buffer.ang[i2] + d2min*buffer.ang[i1] ) / (d1min+d2min)
		vel = ( d1min*buffer.vel[i2] + d2min*buffer.vel[i1] ) / (d1min+d2min)
	end
	
	-- muscle controls
	ang = ang - 3.14159266 
	local u1 = 0.01 - Kp * ang - Kd * vel
	local u2 = 0.01 + Kp * ang + Kd * vel
	muscle1:add_input( limit(u1) )
	muscle2:add_input( limit(u2) )
	
	-- check if it is time to change the external perturbation moment
	if (t - perturb_time) >= (perturb_interval-0.0025) then
		target_body:add_external_moment( 0, 0, -perturb_moment )
		perturb_moment = perturb_amplitude * 2 * (math.random() - 0.5)
		target_body:add_external_moment( 0, 0, perturb_moment )
		perturb_time = t
		-- print a message to the scone messages window
		-- scone.debug( "perturbation moment changed at " .. t )
	end
	
	--local s = string.format('Time:%8.3f  Angle:%7.3f Moment:%7.3f u1:%6.3f u2:%6.3f', t, ang, perturb_moment, u1, u2 )
	--scone.info(s)
	
	-- stop when angle is outside of the range
	if t > 0.1 then
		return (ang > 1) or (ang < -1)
	end
end

function store_data( frame )
	-- store some values for analysis
	frame:set_value( "perturb_moment", perturb_moment )
end

function limit(x)
	if x < 0 then
		return 0
		--	elseif x > 1 then
		--		return 1
	else
		return x
	end
end

