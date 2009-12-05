
function nosewheel_info()
	return {
		lod_pixels={5,50,0},
		bounding_radius = 7,
		materials={'leg','tyre'}
	}
end
function nosewheel_static(lod)
	set_material('leg', .5, .5, .5, 1, .5, .5, .5, 2.0, 0, 0, 0)
	set_material('tyre', .3, .3, .3, 1, 0,0,0, 1, 0, 0, 0)
	use_material('leg')
	local v6 = v(0, 0, 0)
	local v7 = v(0, 3, 0)
	local v8 = v(0, 5, 0)
	local divs = lod*4
	cylinder(divs, v6, v8, v(0,0,1), .4)
	cylinder(divs, v7, v8, v(0,0,1), .5)
	use_material('tyre')
	xref_cylinder(divs, v(.5,5,0), v(1,5,0), v(0,0,1), 1.0)
end

function nosewheelunit_info()
	return {
		bounding_radius = 7,
		materials={'inside', 'matvar2'}
	}
end
function nosewheelunit_static(lod)
	set_material('inside', .2,.2,.2,1, 0,0,0, 1, 0,0,0)
end
function nosewheelunit_dynamic(lod)
	-- flaps
	local v6 = v(1.5, 0, 6)
	local v7 = v(1.5, 0, -1)
	local v8 = v(0, 0, 6)
	local v9 = v(0, 0, -1)
	set_material('matvar2', get_arg_material(2))

	use_material('inside')
	zbias(1, v(0,0,0), v(0,1,0))
	-- flap internal
	xref_quad(v8, v6, v7, v9)
	-- SHould use parameter material(2) here but param materials not done yet
	use_material('matvar2')
	local flap_ang = 0.5*math.pi*math.clamp(3*get_arg(0),0,1)
	local wheel_ang = 0.5*math.pi*math.clamp(1.5*(get_arg(0)-0.34), 0, 1)
	local vrot = 1.5*v(-math.cos(flap_ang), math.sin(flap_ang), 0)
	xref_quad(v7, v6, v6+vrot, v7+vrot)
	xref_quad(v7, v7+vrot, v6+vrot, v6)

	call_model('nosewheel', v(0,0,0), v(1,0,0),
	v(0,math.sin(wheel_ang),math.cos(wheel_ang)), 1.0)
	zbias(0)
end

function mainwheel_info()
	return {
		lod_pixels = {5,50,0},
		bounding_radius = 8,
		materials = {'leg', 'tyre'}
	}
end
function mainwheel_static(lod)
	local v6 = v(0,0,0)
	local v7 = v(0,3,0)
	local v8 = v(0,5,0)
	-- crossbar
	local v13 = v(0, 5, 1.4)
	local v14 = v(0, 5, -1.4)
	local divs = 4*lod
	set_material('leg', .5,.5,.5,1, 1,1,1, 2, 0,0,0)
	use_material('leg')
	cylinder(divs, v6, v8, v(0,0,1), .4)
	cylinder(divs, v7, v8, v(0,0,1), .5)
	cylinder(4, v13, v14, v(1,0,0), .5)
	set_material('tyre', .3,.3,.3,1, 0,0,0, 1, 0,0,0)
	use_material('tyre')
	xref_cylinder(divs, v(.5, 5, 1.1), v(1, 5, 1.1), v(0,0,1), 1)
	xref_cylinder(divs, v(.5, 5, -1.1), v(1, 5, -1.1), v(0,0,1), 1)
end


function mainwheelunit_info()
	return {
		bounding_radius = 7,
		materials={'inside','matvar2'}
	}
end
function mainwheelunit_static(lod)
	set_material('inside', .2,.2,.2,1, 0,0,0, 1, 0,0,0)
end
function mainwheelunit_dynamic(lod)
	-- flaps
	local v6 = v(1.5, 0, 6)
	local v7 = v(1.5, 0, -1)
	local v8 = v(0, 0, 6)
	local v9 = v(0, 0, -1)
	set_material('matvar2', get_arg_material(2))

	use_material('inside')
	zbias(1, v(0,0,0), v(0,1,0))
	-- flap internal
	xref_quad(v8, v6, v7, v9)
	-- SHould use parameter material(2) here but param materials not done yet
	use_material('matvar2')
	local flap_ang = 0.5*math.pi*math.clamp(3*get_arg(0),0,1)
	local wheel_ang = 0.5*math.pi*math.clamp(1.5*(get_arg(0)-0.34), 0, 1)
	local vrot = 1.5*v(-math.cos(flap_ang), math.sin(flap_ang), 0)
	xref_quad(v7, v6, v6+vrot, v7+vrot)
	xref_quad(v7, v7+vrot, v6+vrot, v6)

	call_model('mainwheel', v(0,0,0), v(1,0,0),
	v(0,math.sin(wheel_ang),math.cos(wheel_ang)), 1.0)
	zbias(0)
end

function ladybird_info()
	return {
		lod_pixels = {50,100,200,0},
		bounding_radius = 35,
		materials={'white','engines','matvar0', 'matvar2',
		'engine_inside'}
	}
end

function ladybird_static(lod)

	local v06 = v(-4.0, -5.0, -20.0);
	local v07 = v(4.0, -5.0, -20.0);

	local v08 = v(-6.0, 4.0, -10.0);
	local v09 = v(6.0, 4.0, -10.0);

	local v10 = v(-14.0, -5.0, -10.0);
	local v11 = v(-10.0, 5.0, 10.0);

	local v29 = v(10.0, 5.0, 10.0);
	local v30 = v(14.0, -5.0, -10.0);
	local v31 = v(-10.0, -5.0, 10.0);
	local v32 = v(10.0, -5.0, 10.0);

	local v33 = v(-12.0, 0.0, 10.0);
	local v34 = v(-12.0, 0.0, 13.0);

	--// thruster jets
	local v38 = v(-12.0, 0.0, 13.0);
	local v39 = v(-15.0, -3.0, -9.0);

	local v40 = v(-30.0, -4.0, 9.0);
	local v41 = v(-29.0, -5.5, 9.0);
	local v42 = v(-29.0, -4.0, 9.0);
	local v43 = v(-10.0, 0.0, -11.0);

	xref_thruster(v38, v(0,0,1), 50, true)
	xref_thruster(v43, v(0,0,-1), 25)
	
	set_material('white',.5,.5,.5,1,1,1,1,100)
	use_material('matvar0')
	-- matvar(0)
	quad(v06,v08,v09,v07)
	quad(v09,v08,v11,v29)
	xref_tri(v08,v06,v10)

	local divs = lod*2

	local wingtip_rear = v(30,-5,10)
	local cpoint_rear = v(20,4,10)

	local leadingedge_mid = v(24,-5,-3)
	local tmp = v(5,0,5)
	local cpoint_leadingedge1 = leadingedge_mid - tmp
	local cpoint_leadingedge2 = leadingedge_mid + tmp


	-- body flat side piece
	local normal = ((v29-v09):cross(v30-v09)):norm()
	local cpoint_bodycurve = 0.5*(v29+v30) + 3.0*(v29-v30):cross(normal):norm()
	xref_flat(divs, normal,
		{ v09 },
		{ v29 },
		{ cpoint_bodycurve, v30 }
		)

	-- top wing bulge
	xref_quadric_bezier_quad(divs,divs,
		wingtip_rear, cpoint_leadingedge2, leadingedge_mid,
		cpoint_rear, v(17,5,0), cpoint_leadingedge1,
		v29, cpoint_bodycurve, v30)


	-- rear
	xref_flat(divs, v(0,0,1),
		{ wingtip_rear },
		{ cpoint_rear, v29 },
		{ v32 }
	)
	quad(v29,v11,v31,v32) -- rear
	use_material('matvar2')
	quad(v10,v06,v07,v30)
	quad(v32,v31,v10,v30)
	-- underside of wing
	xref_tri(v30, wingtip_rear, v32)
	-- wing leading edge underside
	xref_flat(divs, v(0,-1,0),
		{ v30 },
		{ cpoint_leadingedge1, leadingedge_mid },
		{ cpoint_leadingedge2, wingtip_rear }
	)

	zbias(1, v33, v(0,0,1))
	set_material('engines',.3,.3,.3,1,.3,.3,.3,20)
	use_material('engines')
	xref_tube(4*lod, v33, v34, v(0,1,0), 2.5, 3.0)
	use_material('engine_inside')
	-- matanim!!
	xref_circle(4*lod, v33, v(0,0,1), v(0,1,0), 2.5)
	-- wheels my friend
	
end
function ladybird_dynamic(lod)
	set_material('matvar0', get_arg_material(0))
	set_material('matvar2', get_arg_material(2))
	set_material('engine_inside', lerp_materials(get_arg(2)*30.0, {0, 0, 0, 1, 0, 0, 0, 10, .5, .5, 1 },
				{0, 0, 0, 1, 0, 0, 0, 10, 0, 0, .5 }))
	if get_arg(0) ~= 0 then
		local v35 = v(0.0, -5.0, -13.0);
		local v36 = v(-15.0, -5.0, 3.0);
		local v37 = v(15.0, -5.0, 3.0);

		zbias(1, v35, v(0,-1,0))
		call_model('nosewheelunit', v35, v(-1,0,0), v(0,-1,0), 1)
		call_model('mainwheelunit', v36, v(-1,0,0), v(0,-1,0), 1)
		call_model('mainwheelunit', v37, v(-1,0,0), v(0,-1,0), 1)
		zbias(0)
	end

end

function __walruswing_info()
	return {
		lod_pixels = { 20, 50, 0 },
		scale = 25,
		bounding_radius = 2.0,
		materials = {'matvar0'}
	}
end
function __walruswing_static(lod)
	-- bottom front
	local v06 = v(0.0, 0.0, 1.0)
	-- bottom back
	local v07 = v(0.0, 0.0, -1.0)
	-- top front
	local v08 = v(0.0, 1.5, 0.0)
	-- top back
	local v09 = v(0.0, 1.5, -1.5)
	use_material('matvar0')
	local bend = v(0.175,0,0)
	xref_quadric_bezier_quad(1,lod*4, v07, 0.5*(v07+v09), v09,
			0.5*(v06+v07)+bend, bend, 0.5*(v08+v09)+bend,
			v06, 0.5*(v06+v08), v08)
	flat(lod*4, v(0,1,0), { 0.5*(v08+v09)+bend, v09 },
			{ 0.5*(v08+v09)-bend, v08 })

end
function __walruswing_dynamic(lod)
	set_material('matvar0', get_arg_material(0))
end
function walrus_info()
	return {
		scale = 1.0,
		bounding_radius = 70,
		materials = {'matvar0', 'text'}
	}
end
function walrus_static(lod)

	local v06 = v(-5.0, 10.0, -30.0)
	-- 6, top four body verts
	local v07 = v(5.0, 10.0, -30.0)
	local v08 = v(-5.0, 10.0, 30.0)
	local v09 = v(5.0, 10.0, 30.0)

	local v10 = v(-11.16025, -0.6698729, -25.0)
	-- 10, right four body verts
	local v11 = v(-6.160254, -9.330127, -35.0)
	local v12 = v(-11.16025, -0.6698729, 35.0)
	local v13 = v(-6.160254, -9.330127, 30.0)

	local v14 = v(11.16025, -0.6698729, -25.0)
	-- 14, left four body verts
	local v15 = v(6.160254, -9.330127, -35.0)
	local v16 = v(11.16025, -0.6698729, 35.0)
	local v17 = v(6.160254, -9.330127, 30.0)

	local v18 = v(-5.0, -0.6698729, -60.0)
	-- 18, front two verts
	local v19 = v(5.0, -0.6698729, -60.0)


	local v20 = v(0.0, 10.0, 0.0)
			-- 20, top wing
	local v21 = v(-1.0, 0.0, 0.0)

	local v22 = v(0.0, 1.0, 0.0)


			-- 23, right wing
	local v24 = v(0.5, -0.8660254, 0.0)
	local v25 = v(-0.8660254, -0.5, 0.0)

			-- 26, left wing
	local v27 = v(0.5, 0.8660254, 0.0)
	local v28 = v(0.8660254, -0.5, 0.0)

	local v29 = v(-0.0, 0.0, 40.0)
				-- 29, main thruster
	local v30 = v(-11.0, 0.0, -35.0)
				-- 30, retro
	local v31 = v(11.0, 0.0, -35.0)

	local v32 = v(-9.0, 5.0, -30.0)
					-- 32, right
	local v33 = v(-12.0, -5.0, 30.0)
	local v34 = v(12.0, -5.0, -30.0)
				-- 34, left
	local v35 = v(9.0, 5.0, 30.0)
	local v36 = v(0.0, 12.0, -30.0)
				-- 36, top
	local v37 = v(0.0, 12.0, 30.0)
	local v38 = v(0.0, -12.0, -30.0)
				-- 38, bottom
	local v39 = v(0.0, -12.0, 30.0)


	local v42 = v(-5.0, 10.0, -30.0)
	-- 6, top four body verts
	local v43 = v(-11.16025, -0.6698729, 35.0)

	use_material('matvar0')
	quad(v07, v06, v08, v09)
	quad(v13, v11, v15, v17)
	xref_quad(v08, v06, v10, v12)
	xref_quad(v12, v10, v11, v13)
	quad(v09, v08, v12, v16)
	quad(v16, v12, v13, v17)
	
	quad(v06, v07, v19, v18)
	quad(v18, v19, v15, v11)
	xref_tri(v06, v18, v10)
	xref_tri(v10, v18, v11)

	thruster(v29, v(0,0,1), 50, true)
	thruster(v30, v(0,0,-1), 35, true)
	thruster(v31, v(0,0,-1), 35, true)
	thruster(v32, v(-1,0,0), 25)
	thruster(v33, v(-1,0,0), 25)
	thruster(v34, v(1,0,0), 25)
	thruster(v35, v(1,0,0), 25)
	thruster(v36, v(0,1,0), 25)
	thruster(v37, v(0,1,0), 25)
	thruster(v38, v(0,-1,0), 25)
	thruster(v39, v(0,-1,0), 25)

	call_model('__walruswing', v20, v(-1,0,0), v(0,1,0), 1.0)
end
function walrus_dynamic(lod)
	local v06 = v(-5.0, 10.0, -30.0)
	local v07 = v(5.0, 10.0, -30.0)
	local v08 = v(-5.0, 10.0, 30.0)
	local v10 = v(-11.16025, -0.6698729, -25.0)
	local v12 = v(-11.16025, -0.6698729, 35.0)
	local v14 = v(11.16025, -0.6698729, -25.0)
	local v16 = v(11.16025, -0.6698729, 35.0)
	local v20 = v(0.0, 10.0, 0.0)
	local v23 = v(-8.660254, -5.0, 0.0)
	local v26 = v(8.660254, -5.0, 0.0)

	local v40 = v(0.0, -9.330127, -30.0)
			-- 40, nosewheel
	local v41 = v(0.0, -9.330127, 13.0)
			-- 41, mainwheel
	local v54 = (v07 - v14):cross(v16 - v14):norm()
	local v55 = (v06 - v08):cross(v12 - v08):norm()

	set_material('matvar0', get_arg_material(0))
	set_material('text', .2,.2,.2,1)
	use_material('text')
	geomflag(0x8000)
	local reg = get_arg_string(0)
	zbias(1, v16, v54)
	text(reg, v16, v54, v(0,0,-1), 10.0, {xoffset=1, yoffset=.3})
	zbias(1, v10, v55)
	text(reg, v10, v55, v(0,0,1), 10.0, {xoffset=.8, yoffset=.3})
	geomflag(0)
	if get_arg(0) > 0 then
		zbias(1, v40, v(0,-1,0))
		call_model('nosewheelunit', v40, v(-1,0,0), v(0,-1,0), 2.0)
		call_model('mainwheelunit', v41, v(-1,0,0), v(0,-1,0), 2.0)
	end
	zbias(0)
	local ang = math.pi - 0.5 + 0.5*get_arg(0)
	local xaxis = v(math.sin(ang), math.cos(ang), 0)
	call_model('__walruswing', v23, xaxis, v(0,0,-1):cross(xaxis), 1.0)
	ang = 0.5 - 0.5*get_arg(0)
	local xaxis = v(math.sin(ang), math.cos(ang), 0)
	call_model('__walruswing', v26, xaxis, v(0,0,-1):cross(xaxis), 1.0)
end

function bigtrader_info()
	return {
		scale=2.0,
		lod_pixels = {25,50,0},
		bounding_radius = 100,
		materials = {'matvar0','gray','text','engine_inside'}
	}
end
function bigtrader_static(lod)
	local v06 = v(4.0, -3.0, -35.0)
--} },			// 6, nose vertices
	local v07 = v(-4.0, -3.0, -35.0)
--} },
	local v08 = v(-1.0, -7.0, -32.0)
--} },		
	local v09 = v(1.0, -7.0, -32.0)
--} },

	local v10 = v(6.0, 8.0, -20.0)
--} },			// 10, nose section back
	local v11 = v(-6.0, 8.0, -20.0)
--} },				// and extrusion area
	local v12 = v(-10.0, 4.0, -20.0)
--} },			
	local v13 = v(-10.0, -4.0, -20.0)
--} },			
	local v14 = v(-6.0, -8.0, -20.0)
--} },
	local v15 = v(6.0, -8.0, -20.0)
--} },
	local v16 = v(10.0, -4.0, -20.0)
--} },			
	local v17 = v(10.0, 4.0, -20.0)
--} },

	-- midpoints
	local v18 = v(0.0, 0.0, -20.0)
--} },			// 18
	local v19 = v(0.0, 0.0, -16.0)
--} },			// 
	local v20 = v(0.0, 0.0, 4.0)
--} },			// 
	local v21 = v(0.0, 0.0, 8.0)
--} },			// 
	local v22 = v(0.0, 0.0, 26.0)
--} },		// 

	local v23 = v(-0.3826834, 0.9238795, 0.0)
--} },		// 23, tube norm

	local v24 = v(-12.5, 2.0, 10.0)
--} },			// 24, top engine
	local v25 = v(-12.5, 2.0, 30.0)
--} },
	local v26 = v(-12.5, 2.0, 13.0)
--} },
	local v27 = v(-12.5, 2.0, 27.0)
--} },

	local v28 = v(-12.0, -5.5, 10.0)
--} },			// 28, bottom engine
	local v29 = v(-12.0, -5.5, 30.0)
--} },
	local v30 = v(-12.0, -5.5, 13.0)
--} },
	local v31 = v(-12.0, -5.5, 27.0)
--} },

	local v32 = v(-10.0, -4.0, -16.0)
--} },			// 32, right text pos
	local v33 = v(10.0, -4.0, 4.0)
--} },			// left text pos

	
	-- 41, extrusion as at 10 but rev order
	local v41 = v(6.0, 8.0, -20.0)
	local v42 = v(-6.0, 8.0, -20.0)
	local v43 = v(-10.0, 4.0, -20.0)
	local v44 = v(-10.0, -4.0, -20.0)
	local v45 = v(-6.0, -8.0, -20.0)
	local v46 = v(6.0, -8.0, -20.0)
	local v47 = v(10.0, -4.0, -20.0)
	local v48 = v(10.0, 4.0, -20.0)

	set_material('text', .20, .20, .20, 1)
	use_material('matvar0')
	quad(v06,v07,v11,v10)
	xref_tri(v07,v12,v11)
	xref_tri(v07,v13,v12)
	xref_tri(v07,v14,v13)
	xref_tri(v07,v08,v14)
	quad(v07,v06,v09,v08)
	quad(v08,v09,v15,v14)
	quad(v10,v11,v14,v15)
	xref_quad(v11,v12,v13,v14)
	extrusion(v19, v20, v(0,1,0), 1.0,
			v41, v42, v43, v44, v45,v46,v47,v48)
	extrusion(v21, v22, v(0,1,0), 1.0,
			v41, v42, v43, v44, v45,v46,v47,v48)
	xref_tube(8, v24, v25, v(0,1,0), 2.0, 2.5)
	xref_tube(8, v28, v29, v(0,1,0), 2.0, 2.5)
	use_material('engine_inside')
	xref_circle(8, v26, v(0,0,-1), v(0,1,0), 2)
	xref_circle(8, v27, v(0,0,1), v(0,1,0), 2)
	xref_circle(8, v30, v(0,0,-1), v(0,1,0), 2)
	xref_circle(8, v31, v(0,0,1), v(0,1,0), 2)
	set_material('gray', .30, .30, .30,1, .10, .10, .10, 10)
	use_material('gray')
	extrusion(v18, v19, v(0,1,0), .85,
			v41, v42, v43, v44, v45,v46,v47,v48)
	extrusion(v20, v21, v(0,1,0), .85,
			v41, v42, v43, v44, v45,v46,v47,v48)

	xref_thruster(v25, v(0,0,1), 30, true)
	xref_thruster(v29, v(0,0,1), 30, true)
	xref_thruster(v24, v(0,0,-1), 20, true)
	xref_thruster(v28, v(0,0,-1), 20, true)
end
function bigtrader_dynamic(lod)
	set_material('matvar0', get_arg_material(0))
	set_material('engine_inside', lerp_materials(get_arg(2)*30.0, {0, 0, 0, 1, 0, 0, 0, 10, .5, .5, 1 },
				{0, 0, 0, 1, 0, 0, 0, 10, 0, 0, .5 }))
	-- 34, gear pos
	local v34 = v(-5.0, -8.0, -13.0)
	local v35 = v(5.0, -8.0, -13.0)
	local v36 = v(-11.5, -8.0, 25.0)
	local v37 = v(11.5, -8.0, 25.0)
	local v38 = v(-11.5, -8.0, 13.0)
	local v39 = v(11.5, -8.0, 13.0)
	-- 40, dish pos
	local v40 = v(-0.05, 8.0, 15.0)
	local leftText = v(-10.0, 0, -6.4)
	local rightText = v(10.0, 0, -6.4)
	local reg = get_arg_string(0)
	-- this means ignore in collision mesh
	geomflag(0x8000)
	use_material('text')
	zbias(1, leftText, v(-1,0,0))
	text(reg, leftText, v(-1,0,0), v(0,0,1), 4, {center=true})
	zbias(1, rightText, v(1,0,0))
	text(reg, rightText, v(1,0,0), v(0,0,-1), 4, {center=true})
	geomflag(0)

	if get_arg(0) > 0 then
		zbias(1, v34, v(0,-1,0))
		call_model('mainwheelunit', v34, v(-1,0,0), v(0,-1,0), .6)
		call_model('mainwheelunit', v35, v(-1,0,0), v(0,-1,0), .6)
		call_model('mainwheelunit', v36, v(-1,0,0), v(0,-1,0), .5)
		call_model('mainwheelunit', v37, v(-1,0,0), v(0,-1,0), .5)
		call_model('mainwheelunit', v38, v(-1,0,0), v(0,-1,0), .5)
		call_model('mainwheelunit', v39, v(-1,0,0), v(0,-1,0), .5)
	end
	--[[
	PTYPE_ZBIAS, 40, 1, 1,
	PTYPE_SUBOBJECT, 0x8000, SUB_DISH, 40, 1, 100, 200,

	PTYPE_ZBIAS, 0x8000, 0, 0,
	--]]
	zbias(0)
end

function interdictor_info()
	return {
		bounding_radius = 100,
		materials = {'matvar0', 'matvar2', 'engine'}
	}
end

function vlerp(t, v1, v2)
	return t*v2 + (1.0-t)*v1
end

function interdictor_static()
	local v06 = v(0.0, 0.0, -35.0)
	--f } },			// 6, nose point
	local v07 = norm(0.0, 1.0, -0.2)
--} },			// nose normal
	local v08 = v(-6.0, 0.0, -18.0)
--} },			// 8, r edge forward mid
	local v09 = norm(-0.2, 1.0, -0.1)
--} },			// norm
	local v10 = v(-12.0, 0.0, 2.0)
--} },		// 10, r edge back mid
	local v11 = norm(-0.2, 1.0, -0.1)
--} },			// norm
	local v12 = v(-7.5, 0.0, 25.0)
--} },		// 12, r edge back
	local v13 = norm(0.0, 1.0, 0.2)
--} },			// norm
	local v14 = v(0.0, 3.0, -15.0)
--} },			// 14, cockpit front
	local v15 = norm(0.0, 1.0, 0.08)
--} },		// norm
	local v16 = v(-1.5, 3.0, -13.5)
--} },			// 16, cockpit right
	local v17 = norm(0.0, 1.0, 0.08)
--} },		// norm
	local v18 = v(0.0, 3.0, -10.0)
--} },			// 18, cockpit back
	local v19 = norm(0.0, 1.0, 0.08)
--} },		// norm
	local v20 = v(1.5, 3.0, -13.5)
--} },		// 20, cockpit left
	local v21 = norm(0.0, 1.0, 0.08)
--} },		// norm

	local v22 = v(-6.0, 3.0, 5.0)
--} },			// 22, inner right
	local v23 = norm(-0.2, 1.0, -0.2)
--} },			// norm
	local v24 = v(0.0, 3.0, 5.0)
--} },			// 24, inner mid
	local v25 = norm(0.2, 1.0, -0.2)
--} },		// norm

	local v26 = v(-2.0, 2.0, -23.0)
--} },			// 26, fwd midpoint
	local v27 = norm(0.0, 1.0, -0.1)
--} },		// norm
	local v28 = v(-5.0, 2.5, -5.0)
--} },			// 28, right midpoint
	local v29 = norm(-0.08, 1.0, -0.04)
--} },		// norm
	local v30 = v(-7.0, 2.0, 14.0)
--} },		// 30, rear right midpoint
	local v31 = norm(-0.04, 1.0, 0.1)
--} },		// norm

	local v32 = v(-3.0, 3.0, -5.0)
--} },			// 32, central midpoint
	local v33 = v(0.0, 4.0, -12.5)
--} },			// 33, cockpit midpoint
	local v34 = v(-3.75, 4.0, 20.0)
--} },		// 34, nacelle midpoint

	local v35 = v(-7.5, 0.0, 30.0)
--} },		// 35, nacelle outer
	local v36 = v(0.0, 0.0, 30.0)
--} },		// 36, nacelle inner

	-- edge tangents
	local v37 = v(6.0, 4.0, 3.0)
--} },		// 37, edge to mid
	local v38 = v(6.0, 0.0, 3.0)
--} },		//
	local v39 = v(0.0, 4.0, -20.0)
--} },			// 39, rear to mid
	local v40 = v(2.5, 0.0, -20.0)
--} },		//

	local v41 = v(0.0, 0.0, -20.0)
--} },			// 41, mid to nose
	local v42 = v(0.0, -4.0, -20.0)
--} },
	local v43 = v(-6.0, 0.0, -3.0)
--} },			// 43, mid to edge
	local v44 = v(-6.0, -4.0, -3.0)
--} },
	local v45 = v(-2.5, 0.0, 20.0)
--} },			// 45, mid to rear
	local v46 = v(0.0, -4.0, 20.0)
--} },

	local v47 = v(-1.5, 0.0, 0.0)
--} },			// 47, cockpit CW tangents
	local v48 = v(1.5, 0.0, 0.0)
--} },
	local v49 = v(0.0, 0.0, -1.5)
--} },
	local v50 = v(0.0, 0.0, 1.5)
--} },
	local v51 = v(0.0, 0.0, -3.5)
--} },			// 51
	local v52 = v(0.0, 0.0, 3.5)
--} },

	local v53 = v(-10.0, 0.0, 20.0)
--} },			// 53, rear edge tangents
	local v54 = v(10.0, 0.0, 0.0)
--} },
	local v55 = v(4.0, 0.0, -10.0)
--} },			// 55, CCW
	local v56 = v(-5.0, 0.0, -10.0)
--} },

	local v57 = v(0.0, 1.5, 0.0)
--} },			// 57, nacelle tangents
	local v58 = v(0.0, -1.5, 0.0)
--} },
	local v59 = v(0.0, 0.0, -12.0)
--} },
	local v60 = v(0.0, 0.0, 12.0)
--} },

	local v61 = v(-3.75, 4.0, 30.0)
--} },			// 61, nacelle rear midpoint
	local v62 = v(-3.0, 0.0, 0.0)
--} },			// and tangents
	local v63 = v(4.0, 0.0, 0.0)
--} },			// 

	-- underside points
	local v64 = v(-5.0, 0.0, -5.0)
--} },			// 64, upper outer vent
	local v65 = v(0.0, 0.0, -5.0)
--} },			// 65, upper inner vent
	local v66 = v(-5.0, -2.0, -3.0)
--} },			// 66, lower outer vent
	local v67 = v(0.0, -2.0, -3.0)
--} },			// 67, lower inner vent
	local v68 = v(-5.0, -2.0, 30.0)
--} },		// 68, nacelle outer underside
	local v69 = v(0.0, -2.0, 30.0)
--} },		// 69, nacelle inner underside
	local v70 = v(-13.0, 0.0, 14.0)
--} },		// 70, rear underside centre
	local v71 = v(-7.5, 0.0, -3.0)
--} },			// 71, vent outer edge 

	local v72 = v(-3.75, 0.7, 30.0)
--} },			// 72, engine midpoint

	local v73 = v(0.0, 0.0, -15.0)
--} },			// 73, nose gear pos
	local v74 = v(-3.75, -2.0, 15.0)
--} },			// 74, rear right gear
	local v75 = v(3.75, -2.0, 15.0)
--} },			// 75, rear left gear

	local v76 = v(-3.75, 0.7, 32.0)
--} },			// 76, engine end

	local v77 = v(-4.5, -0.3, -4.7)
--} },			// 77, retro vent
	local v78 = v(-0.5, -0.3, -4.7)
--} },			// 
	local v79 = v(-4.5, -1.7, -3.3)
--} },			// 
	local v80 = v(-0.5, -1.7, -3.3)
--} },			// 

	-- main & retro thrusters
	local v81 = v(-3.75, 0.7, 32.0)
--} },			// 81
	local v82 = v(3.75, 0.7, 32.0)
--} },			
	local v83 = v(-2.5, -1.0, -5.0)
--} },
	local v84 = v(2.5, -1.0, -5.0)
--} },

	-- vertical thrusters
	local v85 = v(-9.0, 1.5, 10.0)
--} },			// 85
	local v86 = v(-9.0, -0.5, 10.0)
--} },
	local v87 = v(9.0, 1.5, 10.0)
--} },			// 
	local v88 = v(9.0, -0.5, 10.0)
--} },			// 
	local v89 = v(0.0, 3.5, -8.0)
--} },			// 
	local v90 = v(0.0, -0.5, -25.0)
--} },

	-- horizontal thrusters
	local v91 = v(-8.0, 0.0, 28.0)
--} },			// 91
	local v92 = v(8.0, 0.0, 28.0)
--} },
	local v93 = v(-3.5, 0.0, -25.0)
--} },
	local v94 = v(3.5, 0.0, -25.0)
--} },

	-- text norms
	local v95 = norm(-2.0, -2.5, 0.0)
--} },			// 95
	local v96 = norm(2.0, -2.5, 0.0)
--} },
	local v97 = v(5.0, -2.0, 13.5)
--} },		// 97, 98, reg number text points
	local v98 = v(-5.0, -2.0, 13.5)
--} },

	use_material('matvar0')

	local j = v(-2,3.0,-21)
	-- top nose bit
	xref_cubic_bezier_quad(16, 16, v06, vlerp(.25,v06,v08), vlerp(.75,v06,v08), v08,
		vlerp(.25,v06,v14), j, j, vlerp(.25,v08,v16),
		vlerp(.75,v06,v14), j, j, vlerp(.75,v08,v16),
		v14,vlerp(.25,v14,v16),vlerp(.75,v14,v16),v16)
		
	-- side thingies
	j = 0.25*(v10+v22+v16+v08)
	xref_cubic_bezier_quad(16, 16, v08, vlerp(.25,v08,v10), vlerp(.75,v08,v10), v10,
			vlerp(.25,v08,v16), j, j, vlerp(.25,v10,v22),
			vlerp(.75,v08,v16), j, j, vlerp(.75,v10,v22),
			v16, vlerp(.25,v16,v22), vlerp(.75,v16,v22), v22)
	j = 0.333*(v12+v22+v10)
	-- top wings
	xref_cubic_bezier_tri(16, v12, v12+v56, v10-v55, v10,
			vlerp(.25,v12,v22), j, vlerp(.25,v10,v22),
			vlerp(.75,v12,v22), vlerp(.75,v10,v22),
			v22)
	-- 

	--[[
	PTYPE_COMPSMOOTH | RFLAG_XREF, 0x8000, 5, 26, 27, 6, 7,		// front edge
		COMP_HERM_NOTAN, 8, 9,
		COMP_HERMITE, 16, 1, 37, 38,
		COMP_HERMITE, 14, 1, 49, 48,
		COMP_HERMITE, 6, 7, 41, 42,
		COMP_END,
	PTYPE_COMPSMOOTH | RFLAG_XREF, 0x8000, 5, 28, 29, 8, 9,		// mid edge
		COMP_HERM_NOTAN, 10, 11,
		COMP_HERMITE, 22, 1, 37, 38,
		COMP_HERM_NOTAN, 16, 1,
		COMP_HERMITE, 8, 9, 43, 44,
		COMP_END,
	PTYPE_COMPSMOOTH | RFLAG_XREF, 0x8000, 5, 30, 31, 10, 11,		// rear edge
		COMP_HERMITE, 12, 13, 53, 54, 
		COMP_HERMITE, 22, 1, 39, 40,
		COMP_HERMITE, 10, 11, 43, 44,
		COMP_END,
	PTYPE_COMPFLAT | RFLAG_XREF, 0x8000, 5, 32, 1, 16, 1,		// centre
		COMP_HERM_NOTAN, 22, 1,
		COMP_HERMITE, 24, 1, 59, 60,
		COMP_HERM_NOTAN, 18, 1,
		COMP_HERMITE, 16, 1, 47, 51, 
		COMP_END,
	PTYPE_COMPSMOOTH | RFLAG_XREF, 0x8000, 5, 34, 1, 22, 23,		// nacelle
		COMP_HERMITE, 12, 3, 45, 46,
		COMP_HERM_NOTAN, 35, 3,
		COMP_HERMITE, 61, 1, 57, 63,
		COMP_HERMITE, 36, 0, 63, 58,
		COMP_HERM_NOTAN, 24, 25,
		COMP_HERMITE, 22, 23, 59, 60,
		COMP_END,
		--]]
	use_material('matvar2')
	xref_flat(8, v(0,-1,0),
		{ v12+v56, v10-v55, v10 },
		{ v71 }, { v12 })
	--[[
	PTYPE_COMPFLAT | RFLAG_XREF, 0x8000, 5, 70, 4, 12, 4,		// rear underside
		COMP_HERMITE, 10, 4, 56, 55, 
		COMP_LINE, 71, 4,
		COMP_LINE, 12, 4,
		COMP_END,
	--]]
	-- other underside
	xref_quad(v08,v06,v65,v64)
	xref_quad(v08,v64,v71,v10)
	xref_quad(v64,v65,v67,v66)
	xref_tri(v71,v64,v66)
	xref_quad(v71,v66,v68,v12)
	xref_tri(v12,v68,v35)
	xref_quad(v66,v67,v69,v68)

	xref_flat(8, v(0,0,1),
		{ v36+v57, v61-v62, v61 }, 
		{ v61+v62, v35-v58, v35 },
		{ v68 }, { v69 }, { v36 })
--[[
	PTYPE_COMPFLAT | RFLAG_XREF, 7, 5, 72, 2, 36, 2,		// engine back face
		COMP_HERMITE, 61, 2, 57, 62, 
		COMP_HERMITE, 35, 2, 62, 58,
		COMP_LINE, 68, 2,
		COMP_LINE, 69, 2,
		COMP_LINE, 36, 2,
		COMP_END,

	PTYPE_MATFIXED, 30, 30, 30, 30, 30, 30, 200, 0, 0, 0,
	PTYPE_COMPSMOOTH, 0x8000, 5, 33, 1, 16, 0,		// cockpit
		COMP_HERMITE, 18, 2, 52, 48,
		COMP_HERMITE, 20, 0, 48, 51,
		COMP_HERMITE, 14, 5, 49, 47,
		COMP_HERMITE, 16, 3, 47, 50, 
		COMP_END,

		--]]
	zbias(1, v72, v(0,0,1))
	set_material('engine', .30, .30, .30,1, .30, .30, .30, 20)
	use_material('engine')
	xref_tube(12, v72, v76, v(0,1,0), 2.0, 2.5)
	--PTYPE_TUBE | RFLAG_XREF, 8, 12, 72, 76, 1, 250, 200,
	--[[
	PTYPE_MATANIM, AFUNC_THRUSTPULSE,
		0, 0, 0, 0, 0, 0, 100, 50, 50, 100,
		0, 0, 0, 0, 0, 0, 100, 0, 0, 50,
	PTYPE_CIRCLE | RFLAG_XREF, 9, 12, 72, 2, 1, 200,

	PTYPE_ZBIAS, 77, 120, 1,
//	PTYPE_MATFIXED, 30, 30, 30, 0, 0, 0, 100, 0, 0, 0,
	PTYPE_QUADFLAT | RFLAG_XREF, 77, 78, 80, 79,

	PTYPE_MATFIXED, 20, 20, 20, 0, 0, 0, 100, 0, 0, 0,
	PTYPE_ZBIAS, 98, 95, 1,
	PTYPE_TEXT, 0x8000, 0, 98, 95, 2, 0, 200, 250,
	PTYPE_ZBIAS, 97, 96, 1,
	PTYPE_TEXT, 0x8000, 0, 97, 96, 5, 0, 200, 250,

	PTYPE_ZBIAS, 73, 4, 1,
	PTYPE_SUBOBJECT, 0, SUB_NWUNIT, 73, 4, 2, 100,
	PTYPE_ZBIAS, 74, 4, 1,
	PTYPE_SUBOBJECT, 0, SUB_NWUNIT, 74, 4, 2, 64,
	PTYPE_SUBOBJECT, 0, SUB_NWUNIT, 75, 4, 2, 64,
	
	PTYPE_ZBIAS, 0x8000, 0, 0,
	--]]
end
function interdictor_dynamic(lod)
	set_material('matvar0', get_arg_material(0))
	set_material('matvar2', get_arg_material(2))
end

register_models('nosewheel', 'nosewheelunit', 'mainwheel',
'mainwheelunit', 'ladybird', '__walruswing', 'walrus', 'bigtrader',
'interdictor')