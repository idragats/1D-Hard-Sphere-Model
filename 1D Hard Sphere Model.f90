program Hard_Sphere_Model_1D
implicit none
integer,parameter::e=1,Nmc=100000,N=10
real,parameter::s=1
real::Eold,x12old,T,pic,La,Lb,sumEnew,sumEnew2,MEnew,MEnew2,diakimansi,xnew1,x12,Enew,h,x,thesi(N),m,Enew1,Enew2,x1,x2
Logical::Flag
integer::imc,i,j,p,k,L
open  (10,file="ME-T.txt")
open  (20,file="ME-T2.txt")
open  (30,file="ME-T3.txt")
!open  (40,file="ME-x.txt")
open  (50,file="T-Me.txt")
call random_seed()
L=100*s
do T=0.5,5,0.2
  !La=1.2*N*s
  !Lb=5*N*s

!do L=La,Lb

    h=(L)/float(N-1)
     x=-h
    do k=1,N
        x=x+h
	thesi(k)=x
    end do

  sumEnew=0
  sumEnew2=0
  x1=0

do imc=1,2000

  do i=1,N

     call random_number(m)

      x1=(2.0*m-1.0) * 0.1
         Eold = Energeia(i)        
		  xnew1=thesi(i)
          thesi(i) = thesi(i) + x1
             if (thesi(i)>0.and.thesi(i)<L) then
                Enew = Energeia(i)                              
              Call Metrop(Eold,Enew,T,Flag)
                if (Flag==.FALSE.) then
                   thesi(i)=xnew1
				end if
              else
              thesi(i)=xnew1 
			 end if
	 end do
   Enew=0
      do i=1, N
        Enew=enew+Energeia(i)
      end do
 
 Enew=Enew/N/2.0
 sumEnew=sumEnew+Enew
 sumEnew2=sumEnew2+Enew**2
end do

 
sumEnew=0
sumEnew2=0
!x2=0
do imc=1,NMC
  do i=1,N
     call random_number(m)
        x1=(2.0*m-1.0) * 0.1
         Eold =Energeia(i)        
          xnew1=thesi(i)
          thesi(i) = thesi(i) + x1
             if (thesi(i)>0.and.thesi(i)<L) then
                  Enew=Energeia(i)                              
                 Call Metrop(Eold,Enew,T,Flag)
                if (Flag==.FALSE.) then
                    thesi(i)=xnew1
				end if
			 else
               thesi(i)=xnew1 
			 end if
		x1=x1+thesi(i)	 
	end do
!x1=x1/N
!x2=x2+x1
Enew=0
  do i=1,N
     Enew=enew+Energeia(i)
  end do
 
  Enew=Enew/N/2.0
  sumEnew=sumEnew+Enew
  sumEnew2=sumEnew2+Enew**2
end do
print*,T,sumEnew/NMC,(sumEnew2/NMC)-(sumEnew/NMC)**2
write(50,*) T,sumEnew/NMC,(sumEnew2/NMC)-(sumEnew/NMC)**2
 ! print*,real(N)/L, T, sumEnew/NMC, (sumEnew2/NMC)-(sumEnew/NMC)**2
  !if (T==0.5)  then
  !write(10,*) real(N)/L, sumEnew/NMC,(sumEnew2/NMC)-(sumEnew/NMC)**2
  !else if (T==1.5) then
 !write(20,*) real(N)/L, sumEnew/NMC,(sumEnew2/NMC)-(sumEnew/NMC)**2
 !else if (T==5) then  
 !write(30,*) real(N)/L, sumEnew/NMC,(sumEnew2/NMC)-(sumEnew/NMC)**2
 !end if
! write(40,*) x2/NMC,sumEnew/NMC
end do 
close(10)
close(20)
close(30)
close(40)
close(50)
Contains


function Energeia(i) Result(r)
implicit none
real:: r, x12
integer::i,j
 
 r=0.0
do j=1,N
  if (i==j) cycle
    x12 = thesi(j)-thesi(i)
    r=r+4.0*e*((s/x12)**12-(s/x12)**6)

enddo

End function

 

 

 

 

Subroutine Metrop(Eold,Enew,T,Flag)
Implicit none
Real,intent(in)::Eold,Enew,T
Logical,intent(out)::Flag
Real::DE,a
DE=Enew-Eold
if (DE<0.0) then
  Flag=.True.
else
 Call Random_Number(a)
   if (exp(-DE/T)>=a) then
      Flag=.True.
   else
      Flag=.False.
   end if
end if
end subroutine

 

end program