struct NireODESol{uType,tType}
   t::Array{tType,1}
   u::Array{uType,1}
   iter::Array{Float64,1}
   retcode::Symbol
   f
end





function  IRK8(u0, t0, tf, h, n, m, f, p;  initial_interp=true) 
    order=16
    utype = typeof(u0)
    ttype = typeof(t0)
    itermax = 100
    sdt = sign(h)
    (hb, hc, mu, nu) = IRK8Coefficients(h)
    s = length(hb)
    U = Array{utype}(undef, s)
    Uz = Array{utype}(undef, s)
    L = Array{utype}(undef, s)
    L2 = Array{utype}(undef, s)   
    Dmin=Array{utype}(undef,s)
    F = Array{utype}(undef, s)

    for i in 1:s
      Dmin[i] = zeros(eltype(u0), size(u0))
      F[i] = zeros(eltype(u0), size(u0))
    end

    ej=zeros(eltype(u0), size(u0))

    for i in 1:s
        U[i] = zeros(eltype(u0), size(u0))
        Uz[i] = zeros(eltype(u0), size(u0))
        L[i] = zeros(eltype(u0), size(u0))
        L2[i] = zeros(eltype(u0), size(u0))
    end
    uu = Array{typeof(u0)}(undef, n+1)
    tt = Array{ttype}(undef, n+1)

    iters = zeros(n+1)
    uu[1] = copy(u0)
    tt[1] = t0
    tj = [t0, zero(t0)]
    uj = copy(u0)
    j = 0
    cont = true
    while cont
        j += 1
        for i in 1:m
           it = IRK8step!(tj,uj,ej,p,f,h,hb,hc,mu,nu,U,Uz,L,F,Dmin,itermax,initial_interp)
           iters[j+1] += it
        end
        iters[j+1] /= m
        uu[j+1] = uj + ej
        tt[j+1] = tj[1] + tj[2]
        cont = (sdt*tt[j+1] < sdt*tf) && (j<n)
    end
    sol = NireODESol(tt[1:(j+1)],uu[1:(j+1)],iters[1:(j+1)],:Successs,f)
    return sol
  end

function  IRK8(u0, t0, tf, n, m, f, p;  initial_interp=true) 
    order=16
    utype = typeof(u0)
    ttype = typeof(t0)
    itermax = 100
    h = (tf-t0)/(n*m)
    (hb, hc, mu, nu) = IRK8Coefficients(h)
    s = length(hb)
    U = Array{utype}(undef, s)
    Uz = Array{utype}(undef, s)
    L = Array{utype}(undef, s)
    L2 = Array{utype}(undef, s)   
    Dmin=Array{utype}(undef,s)
    F = Array{utype}(undef, s)

    for i in 1:s
      Dmin[i] = zeros(eltype(u0), size(u0))
      F[i] = zeros(eltype(u0), size(u0))
    end

    ej=zeros(eltype(u0), size(u0))

    for i in 1:s
        U[i] = zeros(eltype(u0), size(u0))
        Uz[i] = zeros(eltype(u0), size(u0))
        L[i] = zeros(eltype(u0), size(u0))
        L2[i] = zeros(eltype(u0), size(u0))
    end
    uu = Array{typeof(u0)}(undef, n+1)
    tt = Array{ttype}(undef, n+1)

    iters = zeros(n+1)
    uu[1] = copy(u0)
    tt[1] = t0
    tj = [t0, zero(t0)]
    uj = copy(u0)
    cont = true
    for j in 1:n
        for i in 1:m
           it = IRK8step!(tj,uj,ej,p,f,h,hb,hc,mu,nu,U,Uz,L,F,Dmin,itermax,initial_interp)
           iters[j+1] += it
        end
        iters[j+1] /= m
        uu[j+1] = uj + ej
        tt[j+1] = tj[1] + tj[2]
    end
    sol = NireODESol(tt,uu,iters,:Successs,f)
    return sol
  end





function IRK8step!(ttj,uj,ej,p,f,h,hb,hc,mu,nu,U,Uz,L,F,Dmin,itermax,initial_interp)     
       s = length(hb)
       elems = s*length(uj)
       tj = ttj[1]
       te = ttj[2]
       if initial_interp
          for is in 1:s
            for k in eachindex(uj)
                aux=0.
                for js in 1:s
                    aux+=nu[is,js]*L[js][k]
                end
                U[is][k]=(uj[k]+ej[k])+aux
            end
          end
       else
       for is in 1:s
            @. U[is] = uj + ej
          end
       end  

       for is in 1:s
            f(F[is], U[is], p, tj + hc[is])
            @. L[is] =  hb[is]*F[is]
            Dmin[is] .= Inf
        end


        iter = true # Initialize iter outside the for loop
        plusIt=true

        j=1

        while (j<itermax && iter)

            j+=1
            iter=false
            D0=0

            for is in 1:s

                Uz[is] .= U[is]
                @. U[is] = uj + (ej + mu[is,1]*L[1] + mu[is,2]*L[2] + mu[is,3]*L[3] + mu[is,4]*L[4] +
                         mu[is,5]*L[5] + mu[is,6]*L[6] + mu[is,7]*L[7] + mu[is,8]*L[8])

            end


            for is in 1:s

                eval=false
                for k in eachindex(uj)
                            DY=abs(U[is][k]-Uz[is][k])
                            if DY>0.
                               eval=true
                               if DY< Dmin[is][k]
                                  Dmin[is][k]=DY
                                  iter=true
                               end
                           else
                               D0+=1
                           end
                end


               if eval==true
                  f(F[is], U[is], p,  tj + hc[is])
                  @. L[is] = hb[is]*F[is]
               end
           end


            if (iter==false && D0<elems && plusIt)
                iter=true
                plusIt=false
            else
                plusIt=true
            end

        end # while

        indices = eachindex(uj)

        for k in indices    #Batura konpentsatuaren parekoa
            e0 = ej[k] 
            for is in 1:s
	         e0 += muladd(F[is][k], hb[is], -L[is][k])
            end
            res = Base.TwicePrecision(uj[k], e0)
            for is in 1:s
	       res += L[is][k]
            end
            uj[k] = res.hi
            ej[k] = res.lo
         end  

         
         res = Base.TwicePrecision(tj, te) + h
         ttj[1] = res.hi
         ttj[2] = res.lo
         return  (j)

end




function IRK8Coefficients(h)

"""
Float64

    mu = [0.5 -0.081894963105581419782 0.040042703777945054533 -0.024721345803200289737 0.016976173236371128183 -0.012225914113298097519 0.0087485667691973301174 -0.0054828082532158983753
          1.0818949631055814198 0.5 -0.086958924300832629584 0.044941126302625700184 -0.028759775474749282864 0.020017127636408726943 -0.014074355889166945133 0.0087485667691973301174
          0.95995729622205494547 1.0869589243008326296 0.5 -0.088093838732308249462 0.046165464814799994642 -0.029603873064977914709 0.020017127636408726943 -0.012225914113298097519
          1.0247213458032002897 0.95505887369737429982 1.0880938387323082495 0.5 -0.088347161109827876402 0.046165464814799994642 -0.028759775474749282864 0.016976173236371128183
          0.98302382676362887182 1.0287597754747492829 0.95383453518520000536 1.0883471611098278764 0.5 -0.088093838732308249462 0.044941126302625700184 -0.024721345803200289737
          1.0122259141132980975 0.97998287236359127306 1.0296038730649779147 0.95383453518520000536 1.0880938387323082495 0.5 -0.086958924300832629584 0.040042703777945054533
          0.99125143323080266988 1.0140743558891669451 0.97998287236359127306 1.0287597754747492829 0.95505887369737429982 1.0869589243008326296 0.5 -0.081894963105581419782
          1.0054828082532158984 0.99125143323080266988 1.0122259141132980975 0.98302382676362887182 1.0247213458032002897 0.95995729622205494547 1.0818949631055814198 0.5];

          hc = h*[0.019855071751231884 0.10166676129318664 0.2372337950418355 0.4082826787521751 0.591717321247825 0.7627662049581645 0.8983332387068134 0.9801449282487681]
          hb = h*[0.05061426814518813 0.11119051722668724 0.15685332293894363 0.181341891689181 0.181341891689181 0.15685332293894363 0.11119051722668724 0.05061426814518813]

"""
  
    mu = convert.(typeof(h), [
               # s=1
               parse(BigFloat,"0.5") parse(BigFloat,"-8.1894963105581497136508164735930892e-02") parse(BigFloat,"+4.0042703777945052339969045601553216e-02") parse(BigFloat,"-2.4721345803200374868044581645082913e-02") parse(BigFloat,"+1.6976173236371093026881708761116103e-02 ") parse(BigFloat,"-1.2225914113298206053942531721943548e-02") parse(BigFloat,"+8.7485667691973688117766530139122287e-03") parse(BigFloat,"-5.4828082532158826793409353214950813e-03")
               # s=2
               parse(BigFloat,"+1.0818949631055814971365081647359309e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.6958924300832723329070964616248016e-02") parse(BigFloat,"+4.4941126302625688139830943466131237e-02") parse(BigFloat,"-2.8759775474749310978230557041068536e-02")  parse(BigFloat,"+2.0017127636408709173710417097426705e-02") parse(BigFloat,"-1.4074355889166929145973516652599415e-02") parse(BigFloat,"+8.7485667691973688117766530139122287e-03")
               # s=3
               parse(BigFloat,"+9.5995729622205494766003095439844678e-01") parse(BigFloat,"+1.0869589243008327233290709646162480e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8093838732308313442213871391320316e-02") parse(BigFloat,"+4.6165464814800034116730885592456977e-02")  parse(BigFloat,"-2.9603873064977937463012598212122325e-02") parse(BigFloat,"+2.0017127636408709173710417097426705e-02") parse(BigFloat,"-1.2225914113298206053942531721943548e-02")
               # s=4
               parse(BigFloat,"+1.0247213458032003748680445816450829e+00") parse(BigFloat,"+9.5505887369737431186016905653386876e-01") parse(BigFloat,"+1.0880938387323083134422138713913203e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8347161109827784250707380600804499e-02")  parse(BigFloat,"+4.6165464814800034116730885592456977e-02") parse(BigFloat,"-2.8759775474749310978230557041068536e-02") parse(BigFloat,"+1.6976173236371093026881708761116103e-02")
               # s=5
               parse(BigFloat,"+9.8302382676362890697311829123888390e-01") parse(BigFloat,"+1.0287597754747493109782305570410685e+00 ") parse(BigFloat,"+9.5383453518519996588326911440754302e-01") parse(BigFloat,"+1.0883471611098277842507073806008045e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8093838732308313442213871391320316e-02") parse(BigFloat,"+4.4941126302625688139830943466131237e-02") parse(BigFloat,"-2.4721345803200374868044581645082913e-02")
               # s=6
               parse(BigFloat,"+1.0122259141132982060539425317219435e+00") parse(BigFloat,"+9.7998287236359129082628958290257329e-01")  parse(BigFloat,"+1.0296038730649779374630125982121223e+00")  parse(BigFloat,"+9.5383453518519996588326911440754302e-01") parse(BigFloat,"+1.0880938387323083134422138713913203e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.6958924300832723329070964616248016e-02") parse(BigFloat,"+4.0042703777945052339969045601553216e-02")
               # s=7
               parse(BigFloat,"+9.9125143323080263118822334698608777e-01") parse(BigFloat,"+1.0140743558891669291459735166525994e+00")  parse(BigFloat,"+9.7998287236359129082628958290257329e-01")  parse(BigFloat,"+1.0287597754747493109782305570410685e+00") parse(BigFloat,"+9.5505887369737431186016905653386876e-01 ") parse(BigFloat,"+1.0869589243008327233290709646162480e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.1894963105581497136508164735930892e-02")
               # s=8
               parse(BigFloat,"+1.0054828082532158826793409353214951e+00") parse(BigFloat,"+9.9125143323080263118822334698608777e-01")  parse(BigFloat,"+1.0122259141132982060539425317219435e+00")  parse(BigFloat,"+9.8302382676362890697311829123888390e-01") parse(BigFloat,"+1.0247213458032003748680445816450829e+00")  parse(BigFloat,"+9.5995729622205494766003095439844678e-01") parse(BigFloat,"+1.0818949631055814971365081647359309e+00") parse(BigFloat,"0.5")
        ])

        hb = convert.(typeof(h),
           h*[ parse(BigFloat,"+5.0614268145188129576265677154981094e-02"),
               parse(BigFloat,"+1.1119051722668723527217799721312045e-01"),
               parse(BigFloat,"+1.5685332293894364366898110099330067e-01"),
               parse(BigFloat,"+1.8134189168918099148257522463859781e-01"),
               parse(BigFloat,"+1.8134189168918099148257522463859781e-01"),
               parse(BigFloat,"+1.5685332293894364366898110099330067e-01"),
               parse(BigFloat,"+1.1119051722668723527217799721312045e-01"),
               parse(BigFloat,"+5.0614268145188129576265677154981094e-02")
             ])


        hc= convert.(typeof(h),
           h*[ parse(BigFloat,"+1.9855071751231884158219565715263505e-02"),
               parse(BigFloat,"+1.0166676129318663020422303176208480e-01"),
               parse(BigFloat,"+2.3723379504183550709113047540537686e-01"),
               parse(BigFloat,"+4.0828267875217509753026192881990801e-01"),
               parse(BigFloat,"+5.9171732124782490246973807118009203e-01"),
               parse(BigFloat,"+7.6276620495816449290886952459462321e-01"),
               parse(BigFloat,"+8.9833323870681336979577696823791522e-01"),
               parse(BigFloat,"+9.8014492824876811584178043428473653e-01")
        ])

""" Interpolate coefficients """
        nu = convert.(typeof(h),
           [#s=1
            parse(BigFloat,"-2.3580434359907869272604448668704418e-02") parse(BigFloat,"+3.7564028204350781120220588381324588e-02") parse(BigFloat,"-5.2311682298576367742803308364795570e-02") parse(BigFloat,"+7.2153031436804597582100216023947545e-02") parse(BigFloat,"-1.0367897445499288427629400186099435e-01") parse(BigFloat,"+1.6283634881117387717390478413870362e-01") parse(BigFloat,"-3.0234060551288978571240477869701592e-01") parse(BigFloat,"+7.6796597320813176288773924435846961e-01")
            #s=2
            parse(BigFloat,"-7.3721915347372876383416365681676195e-01") parse(BigFloat,"+1.1683495316711616559345205630057529e+00") parse(BigFloat,"-1.6094743347141669348961276646308222e+00") parse(BigFloat,"+2.1754011196351946282108735964539601e+00") parse(BigFloat,"-3.0074288183636139831253139808157757e+00") parse(BigFloat,"+4.3553757998443780247769815248837769e+00") parse(BigFloat,"-6.6165263753225524986317973726684325e+00") parse(BigFloat,"+9.1860239215521482142367391478777201e+00")
            #s=3
            parse(BigFloat,"-1.2297700288095738198718883246880343e+01") parse(BigFloat,"+1.9333875551079588083957605043596343e+01") parse(BigFloat,"-2.6197717938782395621819055691791508e+01") parse(BigFloat,"+3.4375513292689120907660371301821534e+01") parse(BigFloat,"-4.5083720476972868489541375557064921e+01") parse(BigFloat,"+5.9246620811109433433154163113334023e+01") parse(BigFloat,"-7.4772488278621580963027410652542431e+01") parse(BigFloat,"+7.4720772807309976713169982344387532e+01")
            #s=4
            parse(BigFloat,"-1.2360519192775641529208291530012872e+02") parse(BigFloat,"+1.9271209302089311449325157212657236e+02") parse(BigFloat,"-2.5678715794267672761631251206680453e+02") parse(BigFloat,"+3.2738013789187108730152371652322038e+02") parse(BigFloat,"-4.0945219322597510368323043549108768e+02") parse(BigFloat,"+4.9820118998766368260622671514051444e+02") parse(BigFloat,"-5.5851667346968918682289190242626758e+02") parse(BigFloat,"+4.8118799386789480845212083059902066e+02")
            #s=5
            parse(BigFloat,"-7.6506359873607050490475404918734150e+02") parse(BigFloat,"+1.1843177154904885692635422843492179e+03") parse(BigFloat,"-1.5561007397261703002242559137962282e+03") parse(BigFloat,"+1.9384516434527721535098498400498934e+03") parse(BigFloat,"-2.3391243164654789426290647731129215e+03") parse(BigFloat,"+2.7004970213121135878376040434322658e+03") parse(BigFloat,"-2.8241828008846508731118675847222545e+03") parse(BigFloat,"+2.2683066185998760812106944011486825e+03")
            #s=6
            parse(BigFloat,"-2.9864681431085356613428555972148677e+03") parse(BigFloat,"+4.5981543451008061740598602621277281e+03") parse(BigFloat,"-5.9789738886452262969037360399834883e+03") parse(BigFloat,"+7.3248683186540351736460833184871867e+03") parse(BigFloat,"-8.6243827952223739681429053998929657e+03") parse(BigFloat,"+9.6280058867130293266015441780974974e+03") parse(BigFloat,"-9.6712897035873639047772333270447980e+03") parse(BigFloat,"+7.4939111643640677935209888348628654e+03")
            #s=7
            parse(BigFloat,"-7.5322560258736694963352803649936405e+03") parse(BigFloat,"+1.1555244400412346782160633692641882e+04") parse(BigFloat,"-1.4922008897475610276386022948789031e+04") parse(BigFloat,"+1.8084657457385087410934304600504492e+04") parse(BigFloat,"-2.0967791943855618026873030594152413e+04") parse(BigFloat,"+2.2942168855337640310062689548583393e+04") parse(BigFloat,"-2.2527726562941111873889414492088629e+04") parse(BigFloat,"+1.7129917709645864090891579107793160e+04")
            #s=8
            parse(BigFloat,"-1.2513853497502627647496026877301948e+04") parse(BigFloat,"+1.9160065004504862609943292716902368e+04") parse(BigFloat,"-2.4651395180933566824098540012109580e+04") parse(BigFloat,"+2.9706130282097634507632738110877420e+04") parse(BigFloat,"-3.4168328879509508253154372762217900e+04") parse(BigFloat,"+3.7008535511486677576566037943315160e+04") parse(BigFloat,"-3.5939320855736553497649468561881539e+04") parse(BigFloat,"+2.7086805366868524974942204129999368e+04")
        ])

s = length(hb)    
for i in 1:s
    for j in i+1:s
        mu[i,j] = 1 - mu[j,i]
    end
end    
    
hb1 = (h-sum(hb[2:end-1]))/2
hb[1] = hb1
hb[end] = hb1
    
return (hb, hc, mu, nu)
end
