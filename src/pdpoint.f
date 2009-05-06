
      SUBROUTINE pdpoint(X1,X,Y,ssize,pdx,MEDD,MADD,ANGL11)
C      PARAMETER(ndim=2,NC=2000)
      INTEGER ssize,i,j,k,ssize1,count,k1,k2,ssize2
      INTEGER indx(ssize*(ssize-1)/2),brr(ssize)
      double precision X(ssize),Y(ssize),R(ssize),ACTAN(ssize)
      double precision TEMP(ssize),TEMP2(ssize)
      double precision MED(ssize*(ssize-1)/2)
      double precision MAD(ssize*(ssize-1)/2)
      double precision MEDD(2*ssize*ssize)
      double precision MADD(2*ssize*ssize)
      double precision TEMP1(ssize*(ssize-1)/2)
      double precision ANGL1(ssize*ssize*2)
      double precision ANGL11(ssize*ssize*2)
      double precision MED1(ssize), MAD1(ssize),select
      double precision X1(2),minx,pdx,ax,rx,out,RRABS

      k1=INT((ssize+1)/2)
      k2=INT((ssize+2)/2)
      minx=1
      
      if(x1(1).le.minx) then
         minx=x1(1)
      endif
      do i=1,ssize
         if(X(i).le.minx)then
            minx=X(i)
         endif
         brr(i)=i
      enddo
C
      if(minx.le.0)then
         x1(1)=x1(1)-(minx)+1
         do i=1,ssize
            X(i)=X(i)-(minx)+1
         enddo
      endif
C
      do i=1,ssize
         ACTAN(i)=DATAN(Y(i)/X(i))
         R(i)=DSQRT(X(i)**2+Y(i)**2)
      enddo
      
      ax=DATAN(x1(2)/x1(1))
      rx=DSQRT(x1(2)**2+x1(1)**2)
C        
      out=0
      pdx=0

C     CONSIDER ALL ANGLES FORMED BETWEEN THE POINT X1 AND
C     ALL DATA POINTS

CC    Pointwise Angles
CC
      count=1
      ssize2=ssize-1
      do i=1,ssize2
         do j=i+1,ssize
            if(Y(i).eq.Y(j))then
               ANGL1(count)=3.14159265/2
            else
               ANGL1(count)=DATAN(-(X(i)-X(j))/(Y(i)-Y(j)))
            endif
            count=count+1
         enddo
      enddo
C     Angles formed by X1 and all data points
CC
      do i=1,ssize
          if(Y(i).eq.x1(2))then
              if(X(i).eq.x1(1))then
               ANGL1(count)=0
              else
               ANGL1(count)=3.14159265/2
              endif
          else
             ANGL1(count)=DATAN(-(X(i)-x1(1))/(Y(i)-x1(2)))
          endif
          count=count+1
      enddo                 
CC
      ssize1=ssize*(ssize-1)/2
      do j=1,ssize1
         do i=1,ssize
            TEMP(i)=R(i)*DCOS(ACTAN(i)-ANGL1(j))
            TEMP2(i)=TEMP(i)
         enddo
         MED(j)=select(k1,ssize,temp,brr)/2.0+
     +           select(k2,ssize,temp,brr)/2.0
      
         do i=1,ssize
            TEMP(i)=DABS(TEMP2(i)-MED(j))
            TEMP2(i)=TEMP(i)
         enddo
         MAD(j)=select(k1,ssize,temp,brr)/2.0+
     +           select(k2,ssize,temp,brr)/2.0
         if(MAD(j).eq.0)then
           count=1
           do i=1, ssize2
              do k=i+1,ssize
                 temp1(count)=DABS(temp(i)-temp(k))
                 count=count+1
              enddo
           enddo
           MAD(j)=select(INT((ssize1+1)/2),ssize1,temp1,indx)/2.0
     +           + select(INT((ssize1+2)/2),ssize1,temp1,indx)/2.0
         endif
C         
C     COMPUTER THE OUTLYINGNESS FROM SSIZE1 ANGLES
C
       RRABS=DABS(rx*DCOS(ax-ANGL1(j))-MED(j))/MAD(j)
            if(RRABS.gt.OUT)then
              OUT=RRABS
            endif
      enddo

C
C      ssize1=ssize*(ssize-1)/2
      do j=1,ssize
         do i=1,ssize 
            TEMP(i)=R(i)*DCOS(ACTAN(i)-ANGL1(ssize1+j))
            TEMP2(i)=TEMP(i)
         enddo
C    
         MED1(j)=select(k1,ssize,temp,brr)/2.0+
     +           select(k2,ssize,temp,brr)/2.0  
C
         do i=1,ssize
            TEMP(i)=DABS(TEMP2(i)-MED1(j))
            TEMP2(i)=TEMP(i)
         enddo
C        
         MAD1(j)=select(k1,ssize,temp,brr)/2.0+
     +           select(k2,ssize,temp,brr)/2.0 
         if(MAD1(j).eq.0)then
           count=1
           do i=1, ssize2
              do k=i+1,ssize
                 temp1(count)=DABS(temp(i)-temp(k))
                 count=count+1
              enddo
           enddo        
           MAD1(j)=select(INT((ssize1+1)/2),ssize1,temp1,indx)/2.0
     +            +select(INT((ssize1+2)/2),ssize1,temp1,indx)/2.0     
         endif
C
C     COMPUTER THE OUTLYINGNESS FROM SSIZE1 ANGLES     
C
         RRABS=DABS(rx*DCOS(ax-ANGL1(ssize1+j))-MED1(j))/MAD1(j)
         if(RRABS.gt.out)then
           out=RRABS
         endif
      enddo       

C      pdx=1/(1+out)

CC     Use the medd and madd and angl11 to get out
            
      ssize1=ssize*ssize*2
      do j=1, ssize1
           RRABS=DABS(rx*DCOS(ax-ANGL11(j))-MEDD(j))/MADD(j)
             if(RRABS.gt.OUT)then
              OUT=RRABS
             endif
      enddo
              
      pdx=1/(1+out)

C
      if(minx.le.0)then
         x1(1)=x1(1)+minx-1
         do i=1,ssize
            X(i)=X(i)+minx-1
         enddo
      endif   
C
      RETURN
      END
