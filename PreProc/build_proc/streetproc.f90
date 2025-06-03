   Implicit none
   integer, parameter :: ngc = 26244
   Real,    parameter :: delx = 4./9. *1000.
   integer Gid,Stid,i,ig
   Integer, Dimension(ngc) :: Row,Col, flag
   Real, Dimension(ngc,18) :: fdir
   Real stl,Stl1,stl2,Stdir,Dir1,Dir2
   Real B1, B2


   open(5,file = 'StreetData_input.csv')
   write(6,*)'Gid, Row, Col, Dir1, B1, Dir2, B2'

   fdir = 0.0
   flag=0
  read(5,*)
1   read(5,*,end=69) Gid,Stid,Stdir,StL,Row(gid), Col(gid)
   do i=1,18
     if (Stdir.LT.i*10) then
         fdir(Gid,i) = fdir(gid,i) + stl
         exit
     endif
   enddo
   flag(gid) =1
  go to 1
69  continue

   do ig=1,ngc
     if(flag(ig).eq.0) cycle
     stl1=0.
     dir1=0
     do i=1,18   
       if( fdir(ig,i).gt.stl1) then
         stl1 = fdir(ig,i)
         dir1 = i
       endif
     enddo
     stl2=0.
     dir2=0
     do i=1,18   
       if(i.eq.dir1) cycle
       if( fdir(ig,i).gt.stl2) then
         stL2 = fdir(ig,i)
         dir2 = i
       endif
     enddo
     B1 = delx*delx/stl1
     B2 = delx*delx/stl2
     dir1= dir1*10 - 5
     dir2= dir2*10 - 5
    if(dir2.gt.0) write(6,23)ig,Row(ig),Col(ig), dir1, B1, dir2, B2
23 format(3(I5,','),F6.1,',',F8.1,',',F6.1,',',F8.1)
    enddo

    end

