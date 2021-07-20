#2021/7/18
#Edited by shwh
using Plots
using Fmt
function parameter()
    xmax = Float64(pi)
    xmin = 0.0
    xsp = 100
    dx = (xmax-xmin)/xsp

    ymax = Float64(pi)
    ymin = 0.0
    ysp = 100
    dy = (ymax-ymin)/ysp

    x = collect(xmin:dx:xmax)
    y = collect(ymin:dy:ymax)

    tstep = 600
    dt = 0.01
    c = 1.0
    alpha = (c*dt/dx)
    u = zeros(tstep, length(x), length(y))
    return x,y,dx,dy,dt,alpha,u,tstep 
end

#set up initial condition
function IC(x,y,u,alpha)
    for i = 1:length(x)
        for j = 1:length(y)
            #=
            u[1,i,j] = exp(-10*(i*dx)^2)*exp(-10*(j*dy)^2)*2
            =#
            u[1,i,j] = exp(-10*(i*dx-2)^2)*exp(-10*(j*dy-2)^2)*2
        end
    end
    ubn = []
    push!(ubn,u[1,1,:])
    push!(ubn,u[1,end,:])
    push!(ubn,u[1,:,1])
    push!(ubn,u[1,:,end])

    p1 = surface(x,y,u[1,:,:],color=cgrad(:bwr),title = ("0.0s"),zlims=(-1,1),bg = :lightgray, camera = (30,80),colorbar = false)
    plot(p1)
    savefig("wave_"*string(Int(0))*".png")    
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            u[2,i,j] = u[1,i,j] + 0.5*alpha^2*(u[1,i+1,j] + u[1,i-1,j] + u[1,i,j+1] + u[1,i,j-1] -4*u[1,i,j])
        end
    end
    bound_d(x,y,u,2)    
    #bound_n(x,y,u,2,ubn)    
    return u,ubn
end

#bound
function bound_d(x,y,u,t)
    u[t,1,:] = u[t,2,:]
    u[t,length(x),:] = u[t,length(x)-1,:]
    u[t,:,1] = u[t,:,2]
    u[t,:,length(y)] = u[t,:,length(y)-1]
end

function bound_n(x,y,u,t,ubn)
    u[t,1,:] = ubn[1]
    u[t,length(x),:] = ubn[2]
    u[t,:,1] = ubn[3]
    u[t,:,length(y)] = ubn[4]
end

function main()
    x,y,dx,dy,dt,alpha,u,tstep = parameter()
    u,ubn = IC(x,y,u,alpha)
    #println(ubn)
    for t = 3:tstep
        tt = dt*(t-1)
        for i = 2:length(x)-1
            for j = 2:length(y)-1
                u[t,i,j] = 2*u[t-1,i,j] - u[t-2,i,j] + alpha^2*(u[t-1,i+1,j] + u[t-1,i-1,j] + u[t-1,i,j+1] + u[t-1,i,j-1] -4*u[t-1,i,j])
            end
        end
        bound_d(x,y,u,t) 
        #bound_n(x,y,u,t,ubn)
        if (t-1)%10 == 0
            p2 = surface(x,y,u[t,:,:],color=cgrad(:bwr),title = (((Fmt.format(f"{:.2f}",tt))*"s")),zlims=(-1,1),bg = :lightgray, camera = (30,80),colorbar = false)
            plot(p2)
            savefig("wave_"*string((Int(t-1)))*".png")
        end
    end
end

main()