abstract BansalYaronProblem


function solve_hmj!(byp::BansalYaronProblem ;  
					maxit::Int = 1000000, 
					crit::Float64 = 1e-6, 
					verbose::Bool = true,
					args...)
	if byp.invΔ == 0
		@show mean(byp.V)
		update_value!(byp; args...)
		copy!(byp.V, byp.newV)
		@show mean(byp.V)
	else
		for iter in 1:maxit
			@show iter
			@show mean(byp.V)
			# update newV
			update_value!(byp; args...)
			# check convergence
			distance = chebyshev(vec(byp.newV, vec(byp.V))
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
end


function DataFrame(byp::BansalYaronProblem)
	pd = 1/byp.θ * log(max(byp.V, 0.0)) - log(byp.ρ)
	m = repeat(byp.μs, outer = [length(byp.σs)])
	s2 = repeat(byp.νD^2 * byp.σs, inner = [length(byp.μs)])
	return DataFrame(V = byp.V, pd = pd, m = m, s2 = s2)
end

function loglinearapprox(byp::BansalYaronProblem)
	κ1 = 0.997
	A1 = (1 - 1/byp.ψ)/(1 - κ1 * exp(-byp.κμ))
	A2 = 0.5*((byp.θ - byp.θ / byp.ψ)^2 + (byp.θ * A1 * κ1 * byp.νμ / byp.νD)^2)/(byp.θ * (1-κ1*exp(-byp.κσ)))
	return (A1, A2)
end



function plot(byp::BansalYaronProblem, symbol::Symbol)
	df = DataFrame(byp)
	if symbol == :s2
		condition = find((df[:m] .== byp.μs[1]) | (df[:m] .== byp.μs[50]) | (df[:m] .== byp.μs[end]))
		plot(df[condition, :], x = "s2", y = "pd", color = "m", Geom.line)
	elseif symbol == :m
		condition = find((df[:s2] .== byp.νD^2 * byp.σs[1]) | (df[:s2] .== byp.νD^2 * byp.σs[50]) | (df[:s2] .== byp.νD^2 * byp.σs[end]))
		plot(df[condition, :], x = "m", y = "pd", color = "s2", Geom.line)
	end
end