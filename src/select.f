
C  Numerical recipes in FORTRAN
C  Return the kth smallest value in the array arr(1:n)
C  The input array will be rearranged to have this value
C  in location arr(k), with all smaller elements moved to
C  arr(1:k-1)(in arbitrary order) and all large elements
C  in arr[K+1..n](also in arbitrary order), while making
C  the corresponding rearrangement of the array brr(1:n).
C
      FUNCTION select(k,n,arr,brr)
      INTEGER k,n
      double precision select,arr(n)
      INTEGER brr(n)
      INTEGER i, ir,j,l,mid
      double precision a,b,temp
      l=1
      ir=n
1     if(ir-l.le.1)then
         if(ir-l.eq.1)then
            if(arr(ir).lt.arr(l))then
               temp=arr(l)
               arr(l)=arr(ir)
               arr(ir)=temp
               temp=brr(l)
               brr(l)=brr(ir)
               brr(ir)=temp
            endif
         endif
         select=arr(k)
         return
      else
         mid=(l+ir)/2
         temp=arr(mid)
         arr(mid)=arr(l+1)
         arr(l+1)=temp
         temp=brr(mid)
         brr(mid)=brr(l+1)
         brr(l+1)=temp
         if(arr(l).gt.arr(ir))then
               temp=arr(l)
               arr(l)=arr(ir)
               arr(ir)=temp
               temp=brr(l)   
               brr(l)=brr(ir)
               brr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
               temp=arr(l+1)
               arr(l+1)=arr(ir)
               arr(ir)=temp
               temp=brr(l+1)
               brr(l+1)=brr(ir)
               brr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
               temp=arr(l)
               arr(l)=arr(l+1)
               arr(l+1)=temp
               temp=brr(l)
               brr(l)=brr(l+1)
               brr(l+1)=temp
         endif
         i=l+1
         j=ir
         a=arr(l+1)
         b=brr(l+1)
3        continue
           i=i+1
         if(arr(i).lt.a)goto 3
4        continue
           j=j-1
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         temp=brr(i)
         brr(i)=brr(j)
         brr(j)=temp
         goto 3
5        arr(l+1)=arr(j)
         brr(l+1)=brr(j)
         arr(j)=a
         brr(j)=b
         if(j.ge.k)ir=j-1
         if(j.le.k)l=i
      endif
      goto 1
      END

C  Numerical recipes in FORTRAN 
C  Sort an array arr(1:n) into ascending numerical order 
C  using the Quicksort algorithm, while making the 
C  corresponding rearrangement of the array brr(1:n).  
C  n is input; arr is replaced on output by its sorted
C  rearrangement


      SUBROUTINE quicksort(n,arr,brr)
      INTEGER n,M,NSTACK
      double precision arr(n)
      INTEGER brr(n)
      PARAMETER (M=7,NSTACK=500)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      double precision a,temp
      INTEGER b,temp1
      jstack=0
      l=1
      ir=n
C
C  INSERTION SORT WHEN SUBARRAY SMALL ENOUGH
C
1     if(ir-l.lt.M)then
          do j=l+1,ir
               a=arr(j)
               b=brr(j)
               do i=j-1,l,-1
                    if(arr(i).le.a)goto 2
                    arr(i+1)=arr(i)
                    brr(i+1)=brr(i)
               enddo 
               i=l-1
2              arr(i+1)=a
               brr(i+1)=b
          enddo 
          if(jstack.eq.0)return
C
C  Pop stack and begin a new round of partitioning
C
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
      else
C
C  Choose median of left, center, and right elements as
C  partitioning element a. Also rearrange so that a(l) <=
C  a(l+1) <= a(ir)
C
          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp
          temp1=brr(k)
          brr(k)=brr(l+1)
          brr(l+1)=temp1
          if(arr(l).gt.arr(ir))then
              temp=arr(l)
              arr(l)=arr(ir)
              arr(ir)=temp
              temp1=brr(l)
              brr(l)=brr(ir)
              brr(ir)=temp1     
          endif
          if(arr(l+1).gt.arr(ir))then  
              temp=arr(l+1)
              arr(l+1)=arr(ir)
              arr(ir)=temp
              temp1=brr(l+1)
              brr(l+1)=brr(ir)
              brr(ir)=temp1       
          endif  
          if(arr(l).gt.arr(l+1))then
              temp=arr(l)
              arr(l)=arr(l+1)
              arr(l+1)=temp
              temp1=brr(l)
              brr(l)=brr(l+1)
              brr(l+1)=temp1  
          endif
C          
C  Initialize pointers for partitioning
C
          i=l+1
          j=ir
          a=arr(l+1)
          b=brr(l+1)
C
C  Partitioning element. Beginning of innermost loop.
C  Scan up to find element > a.
C
3         continue
              i=i+1
          if(arr(i).lt.a)goto 3
C
C  Scan down to find element < a.
C
4         continue
              j=j-1
          if(arr(j).gt.a)goto 4  
C
C  Pointers crossed. Exit with partitioning complete. 
C  Exchange elements
C
          if(j.lt.i)goto 5       
          temp=arr(i)
          arr(i)=arr(j)
          arr(j)=temp
          temp1=brr(i)
          brr(i)=brr(j)
          brr(j)=temp1   
          goto 3
C
C  End of innermost loop. Insert partitioning element.
C
5         arr(l+1)=arr(j)
          arr(j)=a
          brr(l+1)=brr(j)
          brr(j)=b       
          jstack=jstack+2
C
C  Push pointers to larger subarray on stack, process smaller
C  subarray immediately.
C
C          if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
          if(ir-i+1.ge.j-1)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
          else
              istack(jstack)=j-1
              istack(jstack-1)=1
              l=i  
          endif
      endif  
      goto 1
C      RETURN
      END
