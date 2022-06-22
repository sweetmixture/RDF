!Helen D Duncan April 2020
module mathsG
contains
  function gaussian(x, b, sigma2) result(gx)
    !Make a gaussian
    double precision, intent(in) :: x, b, sigma2
    double precision :: sigma, gx, prefix, power
    real :: pi
    pi=acos(-1.0); sigma=sqrt(sigma2)
    prefix=1.0/(sigma*sqrt(2.0*pi))
    power=((x-b)/sigma)**2
    power=power*(-0.5)
    prefix=1.0
    gx=prefix*exp(power)
  end function gaussian
end module mathsG

program rdf
  use mathsG
  !A program to get the RDF from a cluster
  !Bin width and broadening (variance) can be edited in this program
  implicit none
  integer :: i, j, natoms, npairs, ipu, ios, cnt, nbins, opu
  double precision, dimension(:), allocatable :: distance
  double precision, dimension(:,:), allocatable :: coords, opdata
  double precision :: binwidth, tmp, upper, lower, x, b, sig2
  double precision :: xyz(3), a1(3), a2(3), dist, g, csum
  character :: atname
  character (len = 255) :: buffer, ipf, name
  logical :: present
  binwidth=0.1!angstrom
  !sig2=0.01   !variance		default
  sig2=0.05   !variance		

  !print *, "Type in the name of the xyz file (if it is bob.xyz then just type bob)"
  read(5,*) name
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ipf=trim(name)//".xyz"
  
  inquire(file=ipf,exist=present)
  if(present) then
     continue
  else
     error stop "Cannot find the input file"
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  open(newunit=ipu,file=ipf)
  ios=0
  read(ipu,'(a)',iostat=ios) buffer
  read(buffer,*) natoms
  allocate(coords(natoms,3))
  read(ipu,'(a)',iostat=ios) buffer
  do i = 1, natoms
     read(ipu,'(a)',iostat=ios) buffer
     if(ios /= 0) then
        error stop "unexpected end of xyz file"
     end if
     read(buffer,*) atname, xyz
     coords(i,:)=xyz
  end do
  close(ipu)

  
  !npairs = natoms!/(2!*([natoms-2]!))
  npairs=0
  do i = 1, natoms
     do j = 1, natoms
        if(i /= j) then
           npairs = npairs+1
        else
        end if
     end do
  end do
  
  allocate(distance(npairs))

  !write(6,'("There are ", i3, " atoms in the xyz file therefore there are ", i9, " different atom...atom distances")') &
  !     natoms, npairs

  
  cnt=0
  do i = 1, natoms
     do j = 1, natoms
        if(i /= j) then
           cnt=cnt+1
           a1(:)=coords(i,:)
           a2(:)=coords(j,:)
           dist=((a1(1)-a2(1))**2)+((a1(2)-a2(2))**2)+((a1(3)-a2(3))**2)
           dist=sqrt(dist)
           distance(cnt)=dist
        end if
     end do
  end do
  
  !number of bins based on the max distance + 1 angstrom
  tmp=ceiling((maxval(distance)) + 1)
  tmp=tmp/binwidth
  nbins=nint(tmp)+1


  do i = 1, npairs
!     write(6,'(i4, f10.4)') i, distance(i)
  end do

  
  allocate(opdata(nbins,5))
  opdata=0
  opdata(1,1)=0.0
  do i = 2, nbins
     opdata(i,1)=opdata((i-1),1)+binwidth
  end do

  ! simple frequency
  do i = 2, nbins-1
     upper=opdata(i+1,1)
     lower=opdata(i,1)
     do j = 1, npairs
        dist=distance(j)
        if(dist >= lower .and. dist < upper) then
           opdata(i,2)=opdata(i,2)+1
        end if
     end do
  end do

  !gaussian
  do i = 1, nbins !x range
     g=0.0; 
     do j = 1, npairs
        tmp=0;
        x=opdata(i,1)
        b=distance(j)
        tmp=gaussian(x,b,sig2)
        g=g+tmp
     end do
     opdata(i,3)=g
  end do
  
  !normalise so area under graph is 1
  csum=0
  do i = 1, nbins
     csum=csum+opdata(i,3)
  end do
  do i = 1, nbins
     opdata(i,4)=opdata(i,3)/csum
  end do

  !normalise by the number of atoms
  do i = 1, nbins
     opdata(i,5)=opdata(i,3)/natoms
  end do
  
  open(newunit=opu,file="output.csv")
  !write(opu,'(a)') "r (angstrom)  ,  RDF (frequency)  ,  smoothed RDF  ,  normalised smoothed RDF , normalised by natoms"
  write(opu,'(a)') "r (angstrom)  ,  normalised by natoms"
  do i = 1, nbins
     write(opu,'(f10.2, " , ", f15.8)') opdata(i,1), opdata(i,5)
  end do
  close(opu)
end program rdf
