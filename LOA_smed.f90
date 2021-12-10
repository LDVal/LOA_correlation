!This is the code that was used to obtain the curves for Fig.4a in the following paper:
!"Effect of degree correlations above the first shell on the percolation transition" Valdez et al. EPL (Europhysics Letters) 96 (3), 38001 (2011)

!Input: 
!N_node: number of nodes
!election: 1 (if the user wants an assortative network) or -1 (if the user wants a disassortative network)
!N_iterat: number of iterations to correlate the network. The Pearson correlation coefficient increases with "N_iterat"
!nrea: number of realizations

!Output:
!"Pinf_...dat": fraction of nodes belonging to the giant compoenent as a function of p for a link percolation process
!"Smed_...dat": mean finite cluster size <s> as a function of p
!"PearsonCorrCoef_...dat": average Pearson correlation coefficient
!"SmedMax_...dat": mean height of the peak of <s> for a network with "N_node" nodes.

module globals
  implicit none
  save
  integer N_node
  integer N_iterat
  integer kmin,kmax
  integer nptos 
  integer gc,max_mass,cluster_number

  real(8) lambda,r
  
  integer,allocatable,dimension(:)::cluster
  integer,allocatable,dimension(:)::Head,Tail    
  integer,allocatable,dimension(:)::edge
  integer,allocatable,dimension(:)::node
  integer,allocatable,dimension(:)::kk
  
  integer(8),allocatable,dimension(:)::mass  
  
  real(8),allocatable,dimension(:)::fraction1
  real(8),allocatable,dimension(:)::Pk,Smed,Pinf


end module globals


module random
	save
	integer::q1=1,q2=104
	integer ir(670)
	real(8)::nmax=2*(2**30-1)+1
	real(8) sig
	integer::sem(670)
end module random




Program Perc
  use globals
  use random

  implicit none
  
  
  integer i,rea,nrea

  integer(8) product
  
  integer election
  real(8) corr
  real(8) corr_media    
  
  character(3) Var1
  character(6) Var2 
  character(6) Var3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Input
  print*,'Number of nodes'
  N_node    =80000!160000,320000,640000
  !read(*,*) N_node
  print*,'Assortative(1),Disassoortative(-1)'
  election  =1
  !read(*,*) election
  print*,'Number of iterations in the LOA algorithm'
  !Increasing "N_iterat" will increase (decrease) the Pearson correlation coefficient for an assortative (disassortative) network 
  N_iterat  =N_node*0.63d0
  !read(*,*) N_iterat
  print*,'Number of realizations'
  nrea      =10000
  !read(*,*) nrea
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  


  nptos =500  
  !parameters for a Poisson degree distribution
  lambda  =2d0
  kmin    =0
  kmax    =20

  allocate(Pk(0:kmax))
  allocate(mass(0:N_node),cluster(N_node),kk(N_node),Node(N_node))
  allocate(fraction1(nptos),Pinf(nptos),Smed(nptos))


  !Degree distribution for an ER network
  Pk=0
  Pk(0)=exp(-lambda)
  product=1
  do i=1,kmax
     product=product*i         
     Pk(i)=1.d0*exp(-lambda)*lambda**(real(i))/real(product)
  enddo
  Pk=Pk/sum(Pk)
 

  call initialize_random

  do i=nptos,1,-1
     fraction1(i)=0.10d0+0.001d0*i
  enddo

 
  Pinf        =0
  corr_media  =0
  Smed        =0d0


  corr_media=0d0
  
  select case(election)
  case(1)
     Var2='assort'
  case(-1)
     Var2='dissor'
  case(0)
     Var2='uncorr'
  end select

  write(Var3,'(I6)') N_node
  call zeros(var3)


  do rea=1,nrea
     print*,rea
     
     open(3,file='NumberRealizations_'//Var2//'_'//'N_'//Var3//'.dat')
        write(3,*) rea
     close(3)
     
     call ConfModel
     
     call rewindLOA(election)     
     !correlation after LOA
     call PearsonCorrelation(corr) 
     print*,"Pearson correlation coefficient",corr 
     corr_media=corr_media+corr
     !Link Percolation
     call LinkPercolation
 
     mass=0
     max_mass=0
     gc=0

     deallocate(edge)

     jj:if(mod(rea,50)==0)then
      open(1,file='Pinf_'//Var2//'_'//'N_'//Var3//'.dat')
      do i=1,nptos
         write(1,*)fraction1(i),Pinf(i)/dble(rea)
      enddo
      close(1)

      open(1,file='Smed_'//Var2//'_'//'N_'//Var3//'.dat')
      
      do i=1,nptos
         write(1,*)fraction1(i),Smed(i)/dble(rea)
      enddo
      close(1)
      
      open(1,file='PearsonCorrCoef_'//Var2//'_'//'N_'//Var3//'.dat')      
      write(1,*)N_node,corr_media/dble(rea)
      close(1)
      
      open(1,file='SmedMax_'//Var2//'_'//'N_'//Var3//'.dat')      
      write(1,*)N_node,maxval(Smed)/rea
      close(1)  
          
     endif jj
     
   enddo

   print*,'PearsonCorrelation ',corr_media/nrea

end Program Perc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  Configuration Model
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine creates three arrays: kk, edge, and node, which contains all the network topology information.
!kk(i) is the degree of node "i".
!The subset edge(node(i):node(i)+kk(i)-1) contains the list of neighbors of node "i".
subroutine ConfModel
  use globals
  use random
  implicit none
  integer i,j,k
  integer m,n,stubmaxpos,stubminpos
  integer counting, nstubs
  integer nstubAux,stubpos1,stubpos2

  real(8) ws

  integer, allocatable, dimension(:)::listStub,kkAux

  real(8),dimension(kmin-1:kmax)::CumulativePk
  
  
  allocate(kkAux(n_node))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution
  
  CumulativePk   =0d0
  ws             =0d0
  do i=kmin,kmax
     ws                = ws+Pk(i)
     CumulativePk(i)   = ws
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Assigning Connectivity to all nodes
2 kk=0 
  do i=1,N_node
     call rand
     do j=kmin,kmax
        if(CumulativePk(j-1)<r.and.r<=CumulativePk(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo

  nstubs=sum(kk)
  
  if(mod(nstubs,2).ne.0) go to 2   !the total number of stubs must be an even number
  
  allocate(listStub(nstubs))  
  allocate(edge(nstubs))

  listStub      =0
  nstubAux      =0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!List of stubs

  do i=1,n_Node
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStub(nstubAux) =i
     enddo
  enddo

  node(1)       =1
  do i=1,n_node-1
     node(i+1)  = node(i)+kk(i)
  enddo

  edge        =0
  kkAux       =0
  counting    =0
  do while(nstubs>0)
3    if(counting==100) then
	deallocate(listStub,edge)
	go to 2
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Randomly choose two stubs
     
     call rand
     stubpos1    =r*nstubs+1
     call rand
     stubpos2    =r*nstubs+1

     m        =listStub(stubpos1)   !the first stub corresponds to node "m"
     n        =listStub(stubpos2)   !the second stub corresponds to node "n"

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Checking if these two stubs can be connected

     if(m.ne.n) then
	do k=1,kkAux(m)
           if(edge(node(m)+k-1)==n)then
              counting     =counting+1
              go to 3   !if nodes "m" and "n" are already connected, we have to choose another pair of stubs
           endif
	enddo
	edge(node(m)+kkAux(m))   =n
	edge(node(n)+kkAux(n))   =m

	kkAux(m)                 =kkAux(m)+1
	kkAux(n)                 =kkAux(n)+1
	stubmaxpos               =max0(stubpos1,stubpos2)
	stubminpos               =min0(stubpos1,stubpos2)
        !updating the list of stubs
	if(stubmaxpos==nstubs) then
           listStub(stubmaxpos)       =listStub(nstubs)
           listStub(stubminpos)       =listStub(nstubs-1)
	else 
           if (stubmaxpos==nstubs-1) then
              listStub(stubminpos)    =listStub(nstubs)
           else
              listStub(stubmaxpos)    =listStub(nstubs)
              listStub(stubminpos)    =listStub(nstubs-1)
           endif
	endif
	nstubs           =nstubs-2
	counting         =0
     else
        counting           =counting+1
        go to 3 !self-loops are forbidden
     endif
  enddo
  deallocate(kkAux,listStub)
end subroutine ConfModel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!Cluster identification
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cluster_id
  use globals
  implicit none
  

  integer mm,i,j,k
  integer secmass
  integer nburn,cs

  integer,dimension(N_node)::ocupp
  integer,dimension(0:N_node+1)::w
  

  cluster       =0
  cluster_number=0
  mass          =0
  max_mass      =0
  secmass       =0
  ocupp         =1


  do i=1,N_node
     !if(kk(i)==0)cycle
     if (cluster(i).eq.0.and.ocupp(i)==1)then
	cluster_number=cluster_number+1
	w(0)=i
	nburn=1
	cs=0
	do while(nburn.ne.0)
           nburn=nburn-1
           j=w(nburn)
           cluster(j)=cluster_number
           cs=cs+1
           do k=0,kk(j)-1
             mm=edge(node(j)+k)
              if(cluster(mm)==0.and.ocupp(mm)==1) then
                 ocupp(mm)=0
                 w(nburn)=mm
                 nburn=nburn+1
              endif
           enddo
	enddo

	Mass(cluster_number)=cs

	if(cs>max_mass)then
           secmass=max_mass
           max_mass=cs
           gc=cluster_number
	endif

	if((cs<max_mass).AND.(secmass<cs)) then
	  secmass=cs
	endif
	
  endif
  
enddo

end subroutine cluster_id



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   PearsonCorrelation
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PearsonCorrelation(corr)
  use globals
  
  implicit none
  integer i
  integer(8) M
  
  real(8) product_kj,medk1,vark1,medk2,vark2
  real(8) corr


  M        =sum(kk)

  allocate(Head(M/2),Tail(M/2))
  call ListofLinks
  product_kj=0d0
  medk1      =0d0
  vark1      =0d0
  medk2      =0d0
  vark2      =0d0
  do i=1,M/2
     product_kj=product_kj+kk(Head(i))*kk(Tail(i))*2
     medk1      =medk1+(kk(Head(i))+kk(Tail(i)))
     medk2      =medk2+(kk(Head(i))+kk(Tail(i)))
     vark1      =vark1+((kk(Head(i)))**2d0+(kk(Tail(i)))**2d0)
     vark2      =vark2+((kk(Head(i)))**2d0+(kk(Tail(i)))**2d0)
  enddo
  product_kj=product_kj/M
  medk1      =medk1/M
  vark1      =vark1/M-medk1**2d0
  medk2      =medk2/M
  vark2      =vark2/M-medk2**2d0  
  
  corr = (product_kj-medk1*medk2)/(vark1*vark2)**0.5d0
  deallocate(Head,Tail)
end subroutine PearsonCorrelation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!      Rewind LOA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rewindLOA(election)
  use globals
  use random

  implicit none
  integer i,j
  integer Aux
  integer v1,v2,v3,v4
  integer pr,pr0,pr1,pr2
  integer nn
  integer sel1,sel2
  integer nod1,nod2,countfail
  integer c13,c14,c23,c24

  integer election

  allocate(Head(sum(kk)/2),Tail(sum(kk)/2))
  call ListofLinks

  nn=sum(kk)/2

  do j=1,N_iterat

1    call rand
     !firstly, we select a link     
     sel1 =r*nn+1
     v1   =Head(sel1)
     v2   =Tail(sel1)

     countfail=0
     !after choosing the link "(v1,v2)", we will try to rewired this link at most 10000 times. 
     !If we fail 10000 times to rewired this link (countfail=10000), we remove it from the list "Head" and "Tail"
2    if(countfail==10000) then
        Head(sel1)=Head(nn)
        Tail(sel1)=Tail(nn)
        nn=nn-1
        countfail=0
        if(nn==2)goto 7
        go to 1
     endif

     !Now we select another link
     call rand
     sel2 =r*nn+1
     v3   =Head(sel2)
     v4   =Tail(sel2)

     if(v4==v1.OR.v4==v2.OR.v3==v1.OR.v3==v2) then         
        countfail=countfail+1
        go to 2
     endif


     c13=0
     c14=0
     c23=0
     c24=0


     pr0= kk(v1)*kk(v2)+kk(v3)*kk(v4)

     !Now we check if node "v1" is already connected to "v4" or "v3". 
     !If v1 is connected to v3, then c13=1. Otherwise, c13=0. Similarly, c14=1 if "v1" is connected to "v4"
     do i=0,kk(v1)-1
        if(Edge(node(v1)+i)==v3) then
           c13=1
        else
           if(Edge(node(v1)+i)==v4) c14=1
        endif
        if(c14+c13==2) exit
     enddo
     !Now we check if node "v2" is already connected to "v4" or "v3". 
     !If v2 is connected to v3, then c23=1. Otherwise, c23=0. Similarly, c24=1 if "v2" is connected to "v4"
     do i=0,kk(v2)-1
        if(Edge(node(v2)+i)==v3) then
           c23=1
        else
           if(Edge(node(v2)+i)==v4) c24=1
        endif
        if(c24+c23==2) exit
     enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !List of cases where rewiring cannot be done because we do not allow multiple links between any pair of nodes

     ! if node "v3" is already connected to "v1" and "v2" then we cannot rewire the links  (c13+c23==2)
     ! if node "v4" is already connected to "v1" and "v2" then we cannot rewire the links  (c14+c24==2)
     ! if node "v1" is already connected to "v3" and "v4" then we cannot rewire the links  (c13+c14==2)
     ! if node "v2" is already connected to "v3" and "v4" then we cannot rewire the links  (c23+c24==2)
     ! the case c13+c14+c23+c24>=3 we cannot rewire the links because it corresponds to the following cases: 
     !i) node "v1" is already connected to "v3" and "v4" and "v2" is also connected to "v3" and/or "v4"
     !ii)node "v2" is already connected to "v3" and "v4" and "v1" is also connected to "v3" and/or "v4"
     if(c13+c23==2.OR.c14+c24==2.OR.c13+c14==2.OR.c23+c24==2.OR.c13+c14+c23+c24>=3) then
        countfail=countfail+1
        go to 2
     endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !List of cases where rewiring can be done


     ! If node "v1" is not connected to nodes "V3" and "v4", and if "v2" is not connected to nodes v3 and v4
     ! then we will choose the configuration that increases the Pearson's coefficient the most, that is,
     ! we choose the configuration that maximazes k_i1*k_j1+k_i2*k_j2
     if(c13+c14+c23+c24==0) then
        pr1=kk(v1)*kk(v3)+kk(v2)*kk(v4)
        pr2=kk(v1)*kk(v4)+kk(v2)*kk(v3)
        if(election==1)then
           pr=max(pr0,pr1,pr2)
           !if pr0>pr1, and pr0>pr2, then we do not rewire
           if(pr==pr0) then
              countfail=countfail+1
              go to 2
           else
              if(pr==pr1) go to 4
              if(pr==pr2) go to 5
           endif
        else
           pr=min(pr0,pr1,pr2)
           if(pr==pr0) then
              countfail=countfail+1
              go to 2
           else
              if(pr==pr1) go to 4
              if(pr==pr2) go to 5
           endif
        endif
     endif

     !If node "v2" is connected to node "v3" but not to "v4" and "v1" is not connected to "v3" nor "v4", 
     !then we will check if reconnecting "v1" to "v3" and "v2" to "v4" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case
     if(c23==1.and.c13+c14+c24==0) then
        select case(election)
        case(1)
           if(kk(v1)*kk(v3)+kk(v2)*kk(v4)>pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v1)*kk(v3)+kk(v2)*kk(v4)<pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if        
        endselect
     endif

     !If node "v2" is connected to node "v4" but not to "v3" and "v1" is not connected to "v3" nor "v4", 
     !then we will check if reconnecting "v1" to "v4" and "v2" to "v3" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case
          
     if(c24==1.and.c13+c14+c23==0) then
        select case(election)
        case(1)
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)>pr0) then
              go to 5
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)<pr0) then
              go to 5
           else
              countfail=countfail+1
              go to 2
           end if
        endselect
     endif     

     !If node "v1" is connected to node "v4" but not to "v3" and "v2" is not connected to "v3" nor "v4", 
     !then we will check if reconnecting "v1" to "v3" and "v2" to "v4" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case
     if(c14==1.and.c13+c23+c24==0) then
        select case(election)
        case(1)
           if(kk(v1)*kk(v3)+kk(v2)*kk(v4)>pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v1)*kk(v3)+kk(v2)*kk(v4)<pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if
        endselect
     endif
     !If node "v1" is connected to node "v3" but not to "v4" and "v2" is not connected to "v3" nor "v4", 
     !then we will check if reconnecting "v1" to "v4" and "v2" to "v3" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case
     if(c13==1.and.c14+c23+c24==0) then
        select case(election)
        case(1)  
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)>pr0) then
              go to 5
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)<pr0) then
              go to 5
           else
              countfail=countfail+1
              go to 2
           end if
        endselect
     endif

     !If node "v1" is connected to node "v3" but not to "v4" and "v2" is connected to "v4" but not "v3",
     !then we will check if reconnecting "v1" to "v4" and "v2" to "v3" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case    				
     if(c24+c13==2.and.c23+c14==0) then
        select case(election)
        case(1)
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)>pr0) then
              go to 5
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v1)*kk(v4)+kk(v2)*kk(v3)<pr0) then
              go to 5
           else
              cycle
           end if
        end select
     endif
     !If node "v1" is connected to node "v4" but not to "v3" and "v2" is connected to "v3" but not "v4",
     !then we will check if reconnecting "v1" to "v3" and "v2" to "v4" increases the Pearson Corr. Coef. when election=1 (assortative case)
     !election=-1 corresponds to the disassortative case    		
     if(c23+c14==2.and.c13+c24==0) then
        select case(election)
        case(1)
           if(kk(v2)*kk(v4)+kk(v3)*kk(v1)>pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if
        case(-1)
           if(kk(v2)*kk(v4)+kk(v3)*kk(v1)<pr0) then
              go to 4
           else
              countfail=countfail+1
              go to 2
           end if
        end select
     endif




4    Tail(sel1)=v3
     Head(sel2)=v2

     do i=0,kk(v1)-1
        if(Edge(node(v1)+i)==v2) then
           Edge(node(v1)+i)=v3
           exit
        endif
     enddo

     do i=0,kk(v2)-1
        if(Edge(node(v2)+i)==v1) then
           Edge(node(v2)+i)=v4
           exit
        endif
     enddo

     do i=0,kk(v3)-1
        if(Edge(node(v3)+i)==v4) then
           Edge(node(v3)+i)=v1
           exit
        endif
     enddo

     do i=0,kk(v4)-1
        if(Edge(node(v4)+i)==v3) then
           Edge(node(v4)+i)=v2
           exit
        endif
     enddo



     go to 6


5    Tail(sel1)=v4
     Tail(sel2)=v2


     do i=0,kk(v1)-1
        if(Edge(node(v1)+i)==v2) then
           Edge(node(v1)+i)=v4
           exit
        endif
     enddo


     do i=0,kk(v2)-1
        if(Edge(node(v2)+i)==v1) then
           Edge(node(v2)+i)=v3
           exit
        endif
     enddo


     do i=0,kk(v3)-1
        if(Edge(node(v3)+i)==v4) then
           Edge(node(v3)+i)=v2
           exit
        endif
     enddo


     do i=0,kk(v4)-1
        if(Edge(node(v4)+i)==v3) then
           Edge(node(v4)+i)=v1
           exit
        endif
     enddo



6    countfail=0

  enddo

7 deallocate(Head,Tail)

end subroutine rewindLOA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! List of Links
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine generates the list of links/edges. 
!Head(i) and Tail(i) are the ends of the i-th link/edge
subroutine ListofLinks
  
  use globals

  implicit none
  integer i,ll,j,k
  integer neig  
  integer, allocatable, dimension(:)::EdgeAux

  ll=0
  allocate(EdgeAux(sum(kk)))
  
  EdgeAux=Edge
  Head=0
  Tail=0
  
  do i=1,N_node   
     if(kk(i)==0)cycle
     do j=0,kk(i)-1
        neig=EdgeAux(node(i)+j)
        if(neig/=0)then
           ll=ll+1
           Head(ll)=neig
           Tail(ll)=i
           do k=0,kk(neig)              
              if(EdgeAux(node(neig)+k)==i) then
                 EdgeAux(node(neig)+k)=0
                 exit 
              endif            
           enddo
        endif
     enddo
  enddo
  deallocate(EdgeAux)
    
end subroutine ListofLinks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!      LinkPercolation
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine LinkPercolation

  use globals
  use random

  implicit none
  integer i,ii,j,k,sel,nn
  integer Nlinks,nod1,nod2,Aux
  integer ptos
  
  real(8) T,S_Aux1,S_Aux2
  
  integer, allocatable, dimension(:)::EdgeAux

  allocate(Head(sum(kk)/2),Tail(sum(kk)/2))
  call ListofLinks
  
  nn    =sum(kk)/2
  Nlinks=nn
  
  do ptos=nptos,1,-1
     T=fraction1(ptos)

     kill:do while(1.d0*nn/dble(Nlinks)>T)
     
        call rand
        sel   =r*nn+1
        nod1  =Head(sel)
        nod2  =Tail(sel)
        Aux   =kk(nod1)
        
        
        do k=0,Aux-1
           ii=edge(node(nod1)+k)
           if(ii==nod2)then
              edge(node(nod1)+k)=edge(node(nod1)+kk(nod1)-1)
              kk(nod1)=kk(nod1)-1
              exit
           endif
        enddo
        
	Aux=kk(nod2)
        do k=0,Aux-1
           ii=edge(node(nod2)+k)
           if(ii==nod1)then
              edge(node(nod2)+k)=edge(node(nod2)+kk(nod2)-1)
              kk(nod2)=kk(nod2)-1
              exit
           endif
        enddo
        
        Head(sel) =Head(nn)
        Tail(sel) =Tail(nn)
        nn=nn-1
        
     enddo kill
     

     call cluster_id
     
     S_Aux1 =0d0
     S_Aux2 =0d0
     do i=1,cluster_number
	if(i==gc) cycle
	S_Aux1=S_Aux1+mass(i)**2
	S_Aux2=S_Aux2+mass(i)
     enddo
 
     if(S_Aux2/=0d0) then 
        Smed(ptos)=Smed(ptos)+S_Aux1/S_Aux2
     endif     
     
     Pinf(ptos)=Pinf(ptos)+real(max_mass)/real(N_node)
       
  enddo
  deallocate(Head,Tail)


end subroutine LinkPercolation



!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)			!llama al reloj
  !hour=11644141
  !print*,hour
  see=hour
  CALL RANDOM_SEED(put=see)		!semilla de la realizacin
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)=r*nmax
     ir(i)=i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random

!****************************************
subroutine rand
  use globals
  use random
  implicit none

1  q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


subroutine zeros(b)
 character*1 b(3)
do i=1,3
   if (b(i).eq.' ') b(i)='0'
end do
end subroutine zeros



