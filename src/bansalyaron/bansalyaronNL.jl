
type BansalYaronProblemNL <: BansalYaronProblem

	# consumption process parameters
	μ::Float64 
	νD::Float64
	κμ::Float64 
	κσ::Float64 
	νμ::Float64 
	νσ::Float64 

	# utility parameters
	ρ::Float64  
	γ::Float64 
	ψ::Float64
	θ::Float64

	# algorithm parameters
	invΔ::Float64         
	invdμ::Float64
	invdσ::Float64 

	# grid    
	μs::Vector{Float64}
	σs::Vector{Float64}

	# storage array
	C::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}
	Ct::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}
	B::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}

	V::Vector{Float64}
	newV::Vector{Float64}    
	u::Vector{Float64}
end



function BansalYaronProblemNL(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5, invΔ = 0.0, μn = 100, σn = 100)

	θ = (1-γ)/(1-1/ψ)

	# create drift grid +-3 sd of stationary distribution
	μn = 100
	μmin = μ - 6 * νμ / sqrt(2*κμ)
	μmax = μ + 6 * νμ / sqrt(2*κμ)
	μs = collect(linspace(μmin, μmax, μn))
	invdμ = (μn - 1)/(μmax - μmin)

	# create volatility +-3 sd otherwise negative volatility
	σn = 100
	σmin = 1.0 - 3 * νσ / sqrt(2*κσ)
	σmax = 1.0 + 3 * νσ / sqrt(2*κσ)
	if σmin < 0
		σmin = 1e-3
		σmax = 2.0
	end
	σs = collect(linspace(σmin, σmax, σn))
	invdσ = (σn - 1)/(σmax - σmin)

	#matrix
	C = spdiagm(
	    (ones(μn*σn), ones(μn*σn-1), ones(μn*σn-1), ones((σn-1)*μn), ones((σn-1)*μn)),
	    (0, 1, -1, μn, - μn)
	)
	Crows = rowvals(C) 
	Cvals = nonzeros(C) 
	fill!(Cvals, zero(Float64))
	ij = zero(Int)
	for σi in 1:σn, μi in 1:μn
		ij += 1
		krange = nzrange(C, ij)
		rows = Crows[krange]
		current =  min(κμ * (μ - μs[μi]) * invdμ, 0.0) - 0.5 * νμ^2 * σs[σi] * invdμ^2
		if μi > 1
			index = searchsortedfirst(rows, ij - 1)
		else
			index = searchsortedfirst(rows, ij)
		end
		Cvals[krange[index]] += current
		
		current =  - max(κμ * (μ - μs[μi]) * invdμ, 0.0) - 0.5 * νμ^2 * σs[σi] * invdμ^2
		if μi < μn
			index = searchsortedfirst(rows, ij + 1)
		else
			index = searchsortedfirst(rows, ij)
		end
		Cvals[krange[index]] += current

		current = min(κσ * (1.0 - σs[σi]) * invdσ, 0.0) - 0.5 * νσ^2 * σs[σi] * invdσ^2
		if σi > 1
			index = searchsortedfirst(rows, ij - μn)
		else
			index = searchsortedfirst(rows, ij)
		end
		Cvals[krange[index]] += current
		
		current = - max(κσ * (1.0 - σs[σi]) * invdσ, 0.0) - 0.5 * νσ^2 * σs[σi] * invdσ^2
		if σi < σn
			index = searchsortedfirst(rows, ij + μn)
		else
			index = searchsortedfirst(rows, ij)
		end
		Cvals[krange[index]] += current

		# so that sum of C by column equal zero
		current = (invΔ
				+ ρ * θ  
				- (1-γ) * μs[μi] 
				+ 0.5 * (1-γ) * γ * νD^2 * σs[σi] 
				+ νμ^2 * σs[σi] * invdμ^2
				+ max(κμ * (μ - μs[μi]) * invdμ, 0.0) - min(κμ * (μ - μs[μi]) * invdμ, 0.0)
				+ νσ^2 * σs[σi] * invdσ^2
				+ max(κσ * (1.0 - σs[σi]) * invdσ, 0.0) - min(κσ * (1.0 - σs[σi]) * invdσ, 0.0)
				)
		index = searchsortedfirst(rows, ij)
		Cvals[krange[index]] += current
	end
	Ct = C'
	B = deepcopy(Ct)
	# initialize value at stationary value
	V = Array(Float64, μn * σn)
	fill!(V, (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ)))
	newV = deepcopy(V)
	u = fill(zero(Float64), μn*σn)
	BansalYaronProblemNL(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, invΔ, invdμ, invdσ, μs, σs, C, Ct, B, V, newV, u)
end


function residual!(byp::BansalYaronProblemNL, x::Vector{Float64}, out::Vector{Float64})
	ij = zero(Int)
	μn = length(byp.μs)
	σn = length(byp.σs)
	Cvals = nonzeros(byp.C)
	Crows = rowvals(byp.C)
	for σi in 1:σn, μi in 1:μn
		ij += 1
		krange = nzrange(byp.C, ij)
		rows = Crows[krange]
		current = zero(Float64)
		for k in krange
			current += Cvals[k] * x[Crows[k]]
		end
		current += - byp.ρ * byp.θ * max(x[ij], 0.0)^(1-1/byp.θ)
		current -= byp.invΔ * byp.V[ij]
		out[ij] = current
	end
end


function gradient!(byp::BansalYaronProblemNL, x::Vector{Float64}, out::Base.SparseMatrix.SparseMatrixCSC{Float64, Int})
	ij = zero(Int)
	μn = length(byp.μs)
	σn = length(byp.σs)
	copy!(nonzeros(byp.B), nonzeros(byp.Ct))
	Bvals = nonzeros(byp.B)
	Brows = rowvals(byp.B)
	ij = zero(Int)
	for σi in 1:σn, μi in 1:μn
		ij += 1
		krange = nzrange(byp.B, ij)
		rows = Brows[krange]
		current = - byp.ρ * byp.θ * (1-1/byp.θ) * max(x[ij], 0.0)^(-1/byp.θ)
		index = searchsortedfirst(rows, ij)
		Bvals[krange[index]] += current
	end
	out.colptr = byp.B.colptr
	out.rowval = byp.B.rowval
	out.nzval = byp.B.nzval
end


function update_value!(byp::BansalYaronProblemNL; method::Symbol = :trust_region, iterations = 10_000)
	μn = length(byp.μs)
	σn = length(byp.σs)
	f! = (x, out) -> residual!(byp, x, out)
	g! = (x, out)-> gradient!(byp, x, out)
	df = DifferentiableSparseMultivariateFunction(f!, g!)
	out = nlsolve(df, byp.V, method = method, show_trace = true, iterations = iterations, xtol = 1e-6, ftol = 1e-15)
	@assert all(x -> x >= zero(Float64), out.zero)
	byp.newV = out.zero
end





