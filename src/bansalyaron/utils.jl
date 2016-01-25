abstract BansalYaronProblem


function solve!(byp::BansalYaronProblem ;  
					maxit::Int = 1000000, 
					crit::Float64 = 1e-6, 
					verbose::Bool = true,
					args...)
	for iter in 1:maxit
		@show iter
		@show mean(byp.V)
		# update newV
		update_value!(byp; args...)
		# check convergence
		distance = chebyshev(vec(byp.newV), vec(byp.V))
		@show distance
		if distance < crit
			if verbose
				println("hmj solved : $(iter) iterations")
			end
			break
		else
			# update V using newV
			(byp.newV, byp.V) = (byp.V, byp.newV)
		end
	end
end


function plot(byp::BansalYaronProblem, symbol::Symbol)
	pd = 1/byp.θ * log(max(byp.V, 0.0)) - log(byp.ρ)
	m = repeat(byp.μs, outer = [length(byp.σs)])
	s2 = repeat(byp.νD^2 * byp.σs, inner = [length(byp.μs)])
	df = DataFrame(V = byp.V, pricedividend = pd, drift = m, volatility = s2)
	if symbol == :s2
		df = df[df[:, :volatility] .< 0.00010, :]
		condition = find((df[:drift] .== byp.μs[1]) | (df[:drift] .== byp.μs[50]) | (df[:drift] .== byp.μs[end]))
		plot(df[condition, :], x = "volatility", y = "pricedividend", color = "drift", Geom.line)
	elseif symbol == :m
		condition = find((df[:volatility] .== byp.νD^2 * byp.σs[1]) | (df[:volatility] .== byp.νD^2 * byp.σs[50]) | (df[:volatility] .== byp.νD^2 * byp.σs[end]))
		plot(df[condition, :], x = "drift", y = "pricedividend", color = "volatility", Geom.line)
	end
end