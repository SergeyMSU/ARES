

module Solver

    contains

        !@cuf attributes(host, device) & 
        subroutine chlld_gd(n_state, al, be, ge, &  ! Для газовой динамики без магнитных полей
                                    w, qqq1, qqq2, &
                                    dsl, dsp, dsc, &
                                    qqq)
            ! n_state = 0-2 - какой метод используем
            ! al,be,ge - нормаль
            ! w - скорость грани
            ! qqq1,qqq2 - переменные с двух сторон
            ! dsl,dsp,dsc - поверхности разрывов
            ! qqq - выходные потоки
            implicit real*8 (a-h,o-z)
            
            real(8), intent(out) :: dsl, dsp, dsc
            real(8), intent(in) :: al, be, ge, w
            integer(4), intent(in) :: n_state
            
            real(8), intent(in) :: qqq1(5), qqq2(5)
            real(8), intent(out) :: qqq(5)
            real(8) :: FR(5),FL(5)
            real(8) ::  FW(5),UL(5),UZ(5),UR(5)
            real(8) ::  UZL(5),UZR(5)
            real(8) ::  UZZL(5),UZZR(5)
            real(8) ::  dq(5)

            real(8) ::  vL(3),vR(3),bL(3),bR(3)
            real(8) ::  vzL(3),vzR(3),bzL(3),bzR(3)
            real(8) ::  vzzL(3),vzzR(3),bzzL(3),bzzR(3)
            real(8) ::  aco(3,3), qv(3), qb(3)
            real(8) ::  x0, x1, x2, x3, x4, x5, x6, x7, x8, x9
            
            integer(4) :: n_disco, i
            
            n_disco = 1
            x0 = 0.0
            x1 = 1.0
            x2 = 2.0
            x3 = 3.0
            x4 = 4.0
            x5 = 5.0
            x6 = 6.0
            x7 = 7.0
            x8 = 8.0
            x9 = 9.0
            
            !dimension qqq(5),qqq1(5),qqq2(5)
            !dimension FR(5),FL(5)
            !dimension FW(5),UL(5),UZ(5),UR(5)
            !dimension UZL(5),UZR(5)
            !dimension UZZL(5),UZZR(5)
            !dimension dq(5)
            !
            !dimension vL(3),vR(3),bL(3),bR(3)
            !dimension vzL(3),vzR(3),bzL(3),bzR(3)
            !dimension vzzL(3),vzzR(3),bzzL(3),bzzR(3)
            !dimension aco(3,3),qv(3),qb(3)

            !data x0,x1,x2,x3,x4,x5,x6,x7,x8,x9/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./
            !data n_disco /1/  ! Выбор нахождения скорости крайних характеристик

            !c-------  n_state=0   - one speed LAX
            !c-------  n_state=1   - two speed LAX (HLL,(Harten-Lax-van-Leer))
            !c-------  n_state=2   - two-state (3 speed) HLLC (Contact Discontinuity)
            !c-------  n_state=3   - multi-state (5 speed) HLLD (All Discontinuity)


            pi=dacos(-x1)
            cpi4=x4*pi
            cpi8=x8*pi
            spi4=dsqrt(cpi4)

            eps=1.D-12
            epsb=1.D-06
            eps_p=1.D-06
            eps_d=1.D-03
            

            ga=x5/x3
            g1=ga-x1


            wv=w


            r1 =qqq1(1)
            u1 =qqq1(2)
            v1 =qqq1(3)
            w1 =qqq1(4)
            p1 =qqq1(5)


            r2 =qqq2(1)
            u2 =qqq2(2)
            v2 =qqq2(3)
            w2 =qqq2(4)
            p2 =qqq2(5)

            ro=(r2+r1)/x2
            au=(u2+u1)/x2
            av=(v2+v1)/x2
            aw=(w2+w1)/x2
            ap=(p2+p1)/x2

            d=0.0
            aco(1,1)=al
            aco(2,1)=be
            aco(3,1)=ge
            
            if(dabs(al).lt.dabs(be).and.dabs(al).lt.dabs(ge))then
            aix=x1
            aiy=x0
            aiz=x0
            elseif(dabs(be).lt.dabs(ge))then
            aix=x0
            aiy=x1
            aiz=x0
            else
            aix=x0
            aiy=x0
            aiz=x1
            endif
            
            aik=aix*al+aiy*be+aiz*ge
            d=dsqrt(x1-aik**2)
            aco(1,2)=(aix-aik*al)/d
            aco(2,2)=(aiy-aik*be)/d
            aco(3,2)=(aiz-aik*ge)/d
            aco(1,3)=(aiy*ge-aiz*be)/d
            aco(2,3)=(aiz*al-aix*ge)/d
            aco(3,3)=(aix*be-aiy*al)/d
            

            do i = 1, 3
            vL(i)=aco(1,i)*u1+aco(2,i)*v1+aco(3,i)*w1
            vR(i)=aco(1,i)*u2+aco(2,i)*v2+aco(3,i)*w2
            enddo

                    cL =dsqrt(ga*p1/r1)
                    
                    cfL=cL
                    ptL=p1

                    
                    cR =dsqrt(ga*p2/r2)
                    cfR=cR
                    ptR=p2

                    cC =dsqrt(ga*ap/ro)
                    cfC=cC
                    vC1=(vL(1)+vR(1))/x2

            if(n_disco.eq.1)then
            SL=min( (vL(1)-cfL),(vC1-cfC) )
            SR=max( (vR(1)+cfR),(vC1+cfC) )
            endif

            !c       SL=min( (vL(1)-cfL),(vR(1)-cfR),(vC1-cfC) )
            !c       SR=max( (vL(1)+cfL),(vR(1)+cfR),(vC1+cfC) )

            if(n_disco.eq.0)then
            SL=min( (vL(1)-cfL),(vR(1)-cfR) )
            SR=max( (vL(1)+cfL),(vR(1)+cfR) )
            endif

            if(n_disco.eq.2)then
            SL_1=min( (vL(1)-cfL),(vC1-cfC) )
            SR_1=max( (vR(1)+cfR),(vC1+cfC) )
            SL_2=min( (vL(1)-cfL),(vR(1)-cfR) )
            SR_2=max( (vL(1)+cfL),(vR(1)+cfR) )
                oo = 0.75_8
                oo1= 1.0_8-oo
            SL= oo*SL_1 + oo1*SL_2
            SR= oo*SR_1 + oo1*SR_2
            endif



                    suR=SR-vR(1)
                    suL=SL-vL(1)
                SM=(suR*r2*vR(1)-ptR+ptL-suL*r1*vL(1)) &
                    /(suR*r2-suL*r1)

                dsl=SL
                dsc=SM
                dsp=SR


            if(n_state.eq.0)then
                TR0=dabs(vL(1)+vR(1))/x2+cfC
                TL0=-TR0
                SR=TR0
                SL=TL0
            endif


            upt1=(u1**2+v1**2+w1**2)/x2

            upt2=(u2**2+v2**2+w2**2)/x2

            e1=p1/g1+r1*upt1
            e2=p2/g1+r2*upt2

            FL(1)=r1*vL(1)
            FL(2)=r1*vL(1)*vL(1)+ptL
            FL(3)=r1*vL(1)*vL(2)
            FL(4)=r1*vL(1)*vL(3)
            FL(5)=(e1+ptL)*vL(1) 

            FR(1)=r2*vR(1)
            FR(2)=r2*vR(1)*vR(1)+ptR
            FR(3)=r2*vR(1)*vR(2)   
            FR(4)=r2*vR(1)*vR(3)  
            FR(5)=(e2+ptR)*vR(1)  

                UL(1)=r1
                UL(5)=e1
                UR(1)=r2
                UR(5)=e2
                do ik=1,3
                UL(ik+1)=r1*vL(ik)
                UR(ik+1)=r2*vR(ik)
                enddo

            do ik=1,5
            UZ(ik)=(SR*UR(ik)-SL*UL(ik)+FL(ik)-FR(ik))/(SR-SL)
            enddo

            !c-------- choise for Bn [=UZ(6)] through fan:
            !       if(id_bn.eq.1)UZ(6)=x0
            
            !c----  
            if(n_state.le.1)then

                    do ik=1,5
                    dq(ik)=UR(ik)-UL(ik)
                    enddo

                    TL=SL
                    TR=SR
                    if(SL.gt.wv)then
                    TL=x0
                    do ik=1, 5
                    FW(ik)=wv*UL(ik)
                    enddo
                    endif
                    if(SL.le.wv.and.wv.le.SR)then
                    do ik=1, 5
                    FW(ik)=wv*UZ(ik)
                    enddo
                    endif
                    if(SR.lt.wv)then
                    TR=x0
                    do ik=1, 5
                    FW(ik)=wv*UR(ik)
                    enddo
                    endif


            a=TR*TL
            b=TR-TL

            qqq(1)=(TR*FL(1)-TL*FR(1)+a*dq(1))/b-FW(1)
            qqq(5)=(TR*FL(5)-TL*FR(5)+a*dq(5))/b-FW(5)
            do ik=2,4
                qv(ik-1)=(TR*FL(ik)-TL*FR(ik)+a*dq(ik))/b-FW(ik)
            enddo

            do i = 1,3
            qqq(i+1)=aco(i,1)*qv(1)+aco(i,2)*qv(2)+aco(i,3)*qv(3)
            enddo

            return
            endif
            !c----
            if(n_state == 2)then
            
            suRm=suR/(SR-SM)
            suLm=suL/(SL-SM)
            rzR=r2*suRm
            rzL=r1*suLm
            vzR(1)=SM
            vzL(1)=SM
            ptzR=ptR+r2*suR*(SM-vR(1))
            ptzL=ptL+r1*suL*(SM-vL(1))
            ptz=(ptzR+ptzL)/x2
            
            vzR(2)=UZ(3)/UZ(1)
            vzR(3)=UZ(4)/UZ(1)
            vzL(2)=vzR(2)
            vzL(3)=vzR(3)
            
            vzR(2)=vR(2)
            vzR(3)=vR(3)
            vzL(2)=vL(2)
            vzL(3)=vL(3)
            
            
            ezR=e2*suRm+(ptz*SM-ptR*vR(1))/(SR-SM)
            ezL=e1*suLm+(ptz*SM-ptL*vL(1))/(SL-SM)
            
            
            vzR(2)=vR(2)
            vzR(3)=vR(3)
            vzL(2)=vL(2)
            vzL(3)=vL(3)
            
                UZL(1)=rzL
                UZL(5)=ezL
                UZR(1)=rzR
                UZR(5)=ezR
                do ik=1,3
                UZL(ik+1)=vzL(ik)*rzL
                UZR(ik+1)=vzR(ik)*rzR
                enddo
            
            if(SL.gt.wv)then
                qqq(1)=FL(1)-wv*UL(1)
                qqq(5)=FL(5)-wv*UL(5)
                do ik=2,4
                    qv(ik-1)=FL(ik)-wv*UL(ik)
                enddo
            endif
            
            if(SL.le.wv.and.SM.ge.wv)then
                qqq(1)=FL(1)+SL*(rzL-r1) -wv*UZL(1)
                qqq(5)=FL(5)+SL*(ezL-e1) -wv*UZL(5)
                do ik=2,4
            qv(ik-1)=FL(ik)+SL*(UZL(ik)-UL(ik))-wv*UZL(ik)
                enddo
            else if(SM .le. wv .and. SR .ge. wv)then
                qqq(1)=FR(1)+SR*(rzR-r2) -wv*UZR(1)
                qqq(5)=FR(5)+SR*(ezR-e2) -wv*UZR(5)
                !         do ik=1,3
                !            qv(ik)= FR(ik + 1) + SR * (UZR(ik + 1) - UR(ik + 1)) - wv * UZR(ik + 1)
                !enddo
                qv(1:3) = FR(2:4) + SR * (UZR(2:4) - UR(2:4)) - wv * UZR(2:4)
            endif
            
            if(SR.lt.wv)then
                qqq(1)=FR(1)-wv*UR(1)
                qqq(5)=FR(5)-wv*UR(5)
                do ik=2,4
                    qv(ik - 1)=FR(ik)-wv*UR(ik)
                enddo
                
                
            endif
            
            
            !do i = 1,3
            !      qqq(i+1)=aco(i,1)*qv(1)+aco(i,2)*qv(2)+aco(i,3)*qv(3)
            !enddo
            
            do i = 2,4
                qqq(i)=aco(i - 1,1)*qv(1)+aco(i - 1,2)*qv(2)+aco(i - 1,3)*qv(3)
                enddo
            
            return
            endif

            return
        end subroutine chlld_gd


end module Solver