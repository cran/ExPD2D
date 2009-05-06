
C  SUBROUTINE PDSAMPLE find the projection depth PDS of sample
C  points w.r.t. a bivariate sample (X,Y) with size SSIZE
C  The left half plane is cutted by NCUT directions
C
      SUBROUTINE pdsample(X,Y,ssize,pds,MED1,MAD1,ANGL11)
      PARAMETER(NC=2000)
      INTEGER ssize,i,j,k,ssize1,count,brr(ssize)
      INTEGER indx((ssize-1)*(ssize-1)) 
      double precision X(ssize),Y(ssize),R(ssize),ACTAN(ssize)
      double precision TEMP(ssize)
      double precision PDS(ssize),OUT(ssize),TEMP2(ssize)
      double precision MED(nc),MAD(nc)
      double precision ANGL1(2*ssize*ssize)
      double precision TEMP1((ssize-1)*(ssize-1))
      double precision ANGL11(2*ssize*ssize)      
      double precision MED1(2*ssize*ssize),MAD1(2*ssize*ssize)
      double precision minx,RRABS,select, ttp(ssize-1)      
      INTEGER k1,k2,ssize2,count1,brr1(ssize-1)
      double precision XX(ssize),YY(ssize)

      k1=INT((ssize+1)/2)
      k2=INT((ssize+2)/2)
      minx=1
      ssize2=ssize-1

C     the following assignment were added            
C      do i=1,ssize-1
c         brr1(i)=i
C      enddo
      
      do i=1,ssize
         if(X(i).le.minx)then
            minx=X(i)
         endif
         brr(i)=i  
      enddo

      if(minx.le.0)then
         do i=1,ssize
            X(i)=X(i)-(minx)+1
         enddo
      endif
                                         
      do i=1,ssize
         ACTAN(i)=DATAN(Y(i)/X(i))
         R(i)=DSQRT(X(i)**2+Y(i)**2)
         OUT(i)=0
         PDS(i)=0
      enddo
C

      count=1
      do i=1,ssize
         count1=1
         do j=1,ssize-1     
C     ssize         
             if(j.ne.i)then
                if(X(i).eq.X(j))then
                   if( (Y(j)-Y(i)).gt.0)then
                     ttp(count1) = 3.14159265/2.0
                   else
                     ttp(count1)=3.0*3.14159265/2.0             
                   endif
                else
                   ttp(count1)=DATAN((Y(j)-Y(i))/(X(j)-X(i)))
                   if( (Y(j).gt.Y(i)).and.(X(j).lt.X(i)) )
     +                  ttp(count1)=ttp(count1)+3.14159265
                   if( (Y(j).lt.Y(i)).and.(X(j).lt.X(i)))
     +                  ttp(count1)=ttp(count1)+3.14159265
                   if( (Y(j).lt.Y(i)).and.(X(j).gt.X(i)))
     +                  ttp(count1)=ttp(count1)+3.14159265*2
                endif
                count1=count1+1
             endif
         enddo
         
         ANGL11(count)=select(INT(ssize/2.0),ssize-1,ttp,brr1)
         if((ANGL11(count).gt.3.1415926/2).and.
     +      (ANGL11(count).le.3.1415926*3/2))then
             ANGL11(count)=ANGL11(count)-3.14159265
         endif
         if((ANGL11(count).gt.3.1415926*3/2).and.
     +      (ANGL11(count).le.3.1415926*2))then
             ANGL11(count)=ANGL11(count)-2*3.14159265
         endif
         
         count=count+1
         ANGL11(count)=select(INT((ssize+1)/2.0),ssize-1,ttp,brr1)
         if((ANGL11(count).gt.3.1415926/2).and.
     +      (ANGL11(count).le.3.1415926*3/2))then
             ANGL11(count)=ANGL11(count)-3.14159265
         endif
         if((ANGL11(count).gt.3.1415926*3/2).and.
     +      (ANGL11(count).le.3.1415926*2))then
             ANGL11(count)=ANGL11(count)-2*3.14159265
         endif         
         count=count+1
      enddo
      
      do i=1,ssize    
        XX(i)=X(i)
        YY(i)=Y(i)
        do j=1,ssize
          if(j.ne.i)then
               XX(j)=X(j)
               YY(j)=Y(j)
               do k=1, ssize
                  if((k.ne.i).and.(k.ne.j))then
                    XX(k)=X(k)
                    YY(k)=Y(k)
                    if( (-(x(k)-x(i))*(y(j)-y(i))
     +                 +(y(k)-y(i))*(x(j)-x(i))).lt.0 )then
                       XX(k)=2*X(i)-X(k)
                       YY(k)=2*Y(i)-Y(k)
                    endif
                  endif     
               enddo
               
            count1=1 
            do k=1,ssize-1
C      ssize
              if(k.ne.j)then
                if(XX(j).eq.XX(k))then
                   if( (YY(k)-YY(j)).gt.0)then
                     ttp(count1)=3.14159265/2
                   else   
                     ttp(count1)=3*3.14159265/2
                   endif
                else
                   ttp(count1)=DATAN((YY(k)-YY(j))/(XX(k)-XX(j)))
                   if( (YY(k).gt.YY(j)).and.(XX(k).lt.XX(j)))
     +                  ttp(count1)=ttp(count1)+3.14159265
                   if( (YY(k).lt.YY(j)).and.(XX(k).lt.XX(j)))
     +                  ttp(count1)=ttp(count1)+3.14159265
                   if( (YY(k).lt.YY(j)).and.(XX(k).gt.XX(j)))
     +                  ttp(count1)=ttp(count1)+3.14159265*2
                 endif
                 count1=count1+1
              endif
            enddo
            
            ANGL11(count)=select(INT(ssize/2.0),ssize-1,ttp,brr1)
            if((ANGL11(count).gt.3.1415926/2).and.
     +        (ANGL11(count).le.3.1415926*3/2))then   
               ANGL11(count)=ANGL11(count)-3.14159265
            endif
            if((ANGL11(count).gt.3.1415926*3/2).and.
     +         (ANGL11(count).le.3.1415926*2))then
               ANGL11(count)=ANGL11(count)-2*3.14159265
            endif
            
            count=count+1
            ANGL11(count)=select(INT((ssize+1)/2.0),ssize-1,ttp,brr1)
            if( (ANGL11(count).gt.3.1415926/2).and.
     +         (ANGL11(count).le.3.1415926*3/2) )then   
               ANGL11(count)=ANGL11(count)-3.14159265
            endif
            if((ANGL11(count).gt.3.1415926*3/2).and.
     +        (ANGL11(count).le.3.1415926*2))then
               ANGL11(count)=ANGL11(count)-2*3.14159265
            endif
            count=count+1
          endif
        enddo
      enddo  
      
       
C      ssize1=ssize*(ssize-1)*2
      ssize1=ssize*ssize*2      
      do j=1,ssize1
         ANGL11(j)=DATAN( -1.0/( DTAN(ANGL11(j)) ) )
         do i=1,ssize
            TEMP(i)=R(i)*DCOS(ACTAN(i)-ANGL11(j))
            TEMP2(i)=TEMP(i)
         enddo
         
         MED1(j)=select(k1,ssize,temp,brr)/2.0+
     +           select(k2,ssize,temp,brr)/2.0

         do i=1,ssize
            TEMP(i)=DABS(TEMP2(i)-MED1(j))
            TEMP2(i)=TEMP(i)
         enddo

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
     +           + select(INT((ssize1+2)/2),ssize1,temp1,indx)/2.0  
         endif

         do i=1,ssize
            RRABS=TEMP2(i)/MAD1(j)
            if(RRABS.gt.OUT(i))then
              OUT(i)=RRABS
            endif
         enddo
      enddo

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

      ssize1=2*ssize*ssize
      
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

         do i=1,ssize
            RRABS=TEMP2(i)/MAD(j)
            if(RRABS.gt.OUT(i))then
              OUT(i)=RRABS
            endif
         enddo
      enddo

      do i=1,ssize
         PDS(i)=1/(1+OUT(i))
      enddo

      if(minx.le.0)then
         do i=1,ssize
            X(i)=X(i)+minx-1
         enddo
      endif
      RETURN
      END


