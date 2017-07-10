c     ***************************
      Subroutine padecof(g,z,p,nn,nom)
c     ***************************
c
c     Recursion for Pade coefficient (J.Serene)
      implicit real*8 (a-h,o-z)
      complex*16 g(nom),z(nom),p(nom,nom)
      do j=1,nn
      p(1,j) = g(j)
      enddo
      do j=2,nn
      do i=2,j
      p(i,j)=(p(i-1,i-1)-p(i-1,j))/(z(j)-z(i-1))/p(i-1,j)
      enddo
      enddo
      return
      end
c     ****************************
      subroutine gpade(ge,e,z,p,nn,nom)
c     ****************************
c
c     Calculation of a Green's function for a given pade-coeff-p(i,j)
c     on the real axis e=e+i0
c     NOTE: added a and b in call statement (for SUN)

      implicit real*8 (a-h,o-z)
      complex*16 ge,e,z(nom),p(nom,nom),a(0:nom),b(0:nom)

      a(0)=0.d0
      a(1)=p(1,1)
      b(0)=1.d0
      b(1)=1.d0
      do i=1,nn-1
      a(i+1)=a(i)+(e-z(i))*p(i+1,i+1)*a(i-1)
      b(i+1)=b(i)+(e-z(i))*p(i+1,i+1)*b(i-1)
      enddo

      ge=a(nn)/b(nn)

      return
      end

