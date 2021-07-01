using FFTW

function FourierCoeffsOdd(f,ta,tb,M)
  N = 2*M+1
  h = (tb-ta)/N  
  tt = range(ta,step=h,length=N)
  ff = f.(tt)
  hf = 1/N*rfft(ff)
  return (tt,ff,hf)
end

function FourierEvalOdd(ta,tb,hf,N)
    M = length(hf)-1
    if (N<2M+1) 
        error("$M=length(hf)-1 should be less than (N+1)/2=$div(N+1,2)")
    end
    h = (tb-ta)/N  
    TT = range(ta,step=h,length=N)
    MM = div(N,2)
    hF = zeros(Complex{Float64},MM+1)
    hF[1:(M+1)] = hf[1:(M+1)]
    PP = N*real.(irfft(hF,N))
    return (TT, PP)
end

function FourierCoeffsEven(f,ta,tb,M)
  N = 2*M
  h = (tb-ta)/N  
  tt = range(ta,step=h,length=N)
  ff = f.(tt)
  hf = 1/N*rfft(ff)
  return (tt,ff,hf)
end

function FourierEvalEven(ta,tb,hf,N)
    M = length(hf)-1
    if (N<2M) 
        error("$M=length(hf)-1 should be less than or equal to N/2=$div(N,2)")
    end
    h = (tb-ta)/N  
    TT = range(ta,step=h,length=N)
    MM = div(N,2)
    hF = zeros(Complex{Float64},MM+1)
    hF[1:M] = hf[1:M]
    hF[M+1] = hf[M+1]/2
    PP = N*real.(irfft(hF,N))
    return (TT, PP)
end

function FourierCoeffs(tt,ff)
    N = length(ff)
    hf = 1/N*rfft(ff)
    return (tt,ff,hf)
end

function FourierEval(ta,tb,hf,N)
    if mod(N, 2) == 0
        return FourierEvalEven(ta,tb,hf,N)
    else 
        return FourierEvalOdd(ta,tb,hf,N)
    end
end

function DiffMatrix(N)
    D = zeros(N,N)
    isodd(N) ? auxfcn = csc : auxfcn = cot
    for j in 1:N
        for k in 1:j-1
            D[j,k] = 0.5*(-1)^(j-k)*auxfcn((j-k)*pi/N)
            D[k,j] = -D[j,k]
        end
    end
    return D
end