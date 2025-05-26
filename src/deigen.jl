export deigen
export order

struct DEigenSystem{T,N}
  values_list::Matrix{T}
  vectors_list::Array{T,3}
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, s::DEigenSystem{T,N}) where {T,N}
  summary(io,s); println(io)
  println(io, "values_list:")
  show(io, mime, s.values_list)
  println(io, "\n\nvectors_list:")
  show(io, mime, s.vectors_list)
end

# iteration for upacking into a tuple of variables
Base.iterate(s::DEigenSystem) = (s.values_list, Val(:vectors_list))
Base.iterate(s::DEigenSystem, ::Val{:vectors_list}) = (s.vectors_list, Val(:done))
Base.iterate(s::DEigenSystem, ::Val{:done}) = nothing

order(s::DEigenSystem{T,N}) where {T,N} = N

function (s::DEigenSystem)(inds)
  get_order(s, inds)
end

function get_order(s::DEigenSystem, inds)
  if all(0 .<= inds .<= order(s))
    return _get_order_impl(s, inds)
  else
    error("Index $inds for accessing the orders is out of range. Valid range: 0 to $(order(s))")
  end
end

_get_order_impl(s::DEigen, i::Int) = (s.values_list[:,i+1], s.vectors_list[:,:,i+1])

function _get_order_impl(s::DEigenSystem, inds::AbstractVector{<:Integer})
  shifted_inds = inds .+ 1
  return (s.values_list[:,shifted_inds], s.vectors_list[:,:,shifted_inds])
end

function _get_order_impl(s::DEigenSystem, r::StepRange)
  shifted_r = (first(r)+1):step(r):(stop(r)+1)
  return (s.values_list[:,shifted_r], s.vectors_list[:,:,shifted_r])
end

function deigen(A_list::Vararg{Matrix{T},N}) where {T,N}
  vals0, vecs0 = eigen(A_list[1])
  return deigen(vals0, vecs0, A_list)
end

function deigen(vals0::Vector{T}, vecs0::Matrix{T}, A_list::Vararg{Matrix{T},N}) where {T,N}
  return deigen(vals0, vecs0, A_list)
end

# method for single eigenvalue and single eigenvector
function deigen(val0::T, vec0::Vector{T}, A_list::Vararg{Matrix{T},N}) where {T,N}
  return deigen([val0], reshape(vec0, :, 1), A_list)
end
  
function deigen(vals0::Vector{T}, vecs0::Matrix{T}, A_list::NTuple{N,Matrix{T}}) where {T,N}
  m, n = size(vecs0)
  vals_list = zeros(n,N)
  vecs_list = zeros(m,n,N)
  vals_list[:,:,1] .= vals0
  vecs_list[:,:,1] = vecs0
  b = zeros(m+1)
  E = zeros(m+1,m+1)
  b_hd = zero(T)
  b_tl = zeros(T,m)
  b = zeros(T, m+1)
  x = zeros(T, m+1)
  for i in 1:n
    v0 = view(vecs0,:,i)
    E[1,1] = 0
    E[1,2:end] .= v0
    E[2:end,1] .= v0
    E[2:end,2:end] = vals0[i]*I - A_list[1]
    F = lu!(E)
    for k in 1:(N-1)
      b_hd = 0
      b_tl .= 0
      for l in 0:(k-1)
        b_tl = b_tl + binomial(k,l)*(A_list[k-l+1]*vecs_list[:,i,l+1])
        if l > 1
          b_tl = b_tl - binomial(k,l)*vecs_list[:,i,k-l+1]*vals_list[i,l+1]
          if l < k - 1
            b_hd = b_hd + binomial(k-1,l-1)*vecs_list[:,i,k-l+1]'*vecs_list[:,i,l+1]
          end
        end
      end
      b[1] = -b_hd/2
      b[2:end] = b_tl 
      ldiv!(x,F,b)
      vals_list[i,k+1] = x[1]
      vecs_list[:,i,k+1] .= x[2:end]
    end
  end
  return DEigenSystem{T,N-1}(vals_list, vecs_list)
end
