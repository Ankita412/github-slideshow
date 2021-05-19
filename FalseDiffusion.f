       program FalseDiffusion
       implicit none
       integer:: i, j, l
       INTEGER, PARAMETER:: Nx=102, Ny=102
       REAL(KIND=8):: T(Nx,Ny), aW, aS, aP, Ti(Nx,Ny), dx,dy,tol(Nx*Ny)
     & ,err, x(Nx), y(Ny), aWW, aSS, h, aN, aE
       REAL, PARAMETER:: V=2.0, Length=10.0, Breadth=10.0, dt = 0.0001

       dx=Length/(Nx-2)
       dy=Breadth/(Ny-2)


       do i=1,Nx
       T(i,1)=0.0
       enddo

       do j=1,Ny
       T(1,j)=100.0
       enddo
       h = V*dt/dx
       do i=2,Nx-1
       do j=2,Ny-1
       T(i,j)=0.0
       enddo
       enddo

       tol = 0.0
       do
       Ti = T
       do i=2,Nx-1
       do j=2,Ny-1

       if (((i==2).and.(j>=4)).and.(j/=Ny-1))then
       aW=(4.0/3.0)*h
       aE=(-1.0/3.0)*h
       aS=(7.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=0.0
       aSS=(-1.0/8.0)*h
       elseif (((i==3).and.(j>=4)).and.(j/=Ny-1))then
       aW=(9.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(7.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=(-1.0/3.0)*h
       aSS=(-1.0/8.0)*h
       elseif (((i>=4).and.(i/=Nx-1)).and.(j==2))then
       aW=(7.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(4.0/3.0)*h
       aN=(-1.0/3.0)*h
       aWW=(-1.0/8.0)*h
       aSS=0.0
       elseif (((i>=4).and.(i/=Nx-1)).and.(j==3))then
       aW=(7.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(9.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=(-1.0/8.0)*h
       aSS=(-1.0/3.0)*h
       elseif ((i==2).and.(j==2))then
       aW=(4.0/3.0)*h
       aE=(-1.0/3.0)*h
       aS=(4.0/3.0)*h
       aN=(-1.0/3.0)*h
       aWW=0.0
       aSS=0.0
       elseif ((i==3).and.(j==2))then
       aW=(9.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(4.0/3.0)*h
       aN=(-1.0/3.0)*h
       aWW=(-1.0/3.0)*h
       aSS=0.0
       elseif ((i==2).and.(j==3))then
       aW=(4.0/3.0)*h
       aE=(-1.0/3.0)*h
       aS=(9.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=0.0
       aSS=(-1.0/3.0)*h
       elseif ((i==3).and.(j==3))then
       aW=(9.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(9.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=(-1.0/3.0)*h
       aSS=(-1.0/3.0)*h
       elseif (((i==Nx-1).and.(j>=2)).and.(j/=Ny-1))then
       aW=h
       aE=0.0
       aS=h
       aN=0.0
       aWW=0.0
       aSS=0.0
       elseif ((i>=2).and.(j==Ny-1))then
       aW=h
       aE=0.0
       aS=h
       aN=0.0
       aWW=0.0
       aSS=0.0
       else
       aW=(7.0/8.0)*h
       aE=(-3.0/8.0)*h
       aS=(7.0/8.0)*h
       aN=(-3.0/8.0)*h
       aWW=(-1.0/8.0)*h
       aSS=(-1.0/8.0)*h
       endif
       aP = (1-(aW+aS+aE+aN+aWW+aSS))

       T(i,j) = aW*Ti(i-1,j)+aS*Ti(i,j-1)+aE*Ti(i+1,j)+aN*Ti(i,j+1)+
     & aWW*Ti(i-2,j)+aSS*Ti(i,j-2)+aP*Ti(i,j)
       T(i,Ny)=T(i,Ny-1)
       T(Nx,j)=T(Nx-1,j)
       enddo
       enddo
       T(Nx,j)=T(Nx-1,j)
       
       l=1
       do i=2,Nx-1
       do j=2,Ny-1
       tol(l)=abs(T(i,j)-Ti(i,j))
       l=l+1
       enddo
       enddo
       err=maxval(tol)
       if (err<=1.0E-06) exit
       enddo
       do i=1,Nx
       do j=1,Ny
       print*, "i=",i, "j=",j, "T=",T(i,j)
       enddo
       enddo
       pause
       end

