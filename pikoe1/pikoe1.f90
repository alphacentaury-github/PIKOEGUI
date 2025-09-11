!****************************************************************************************
!*****                                                                              *****
!*****                                 program pikoe                                *****
!*****                                                                              *****
!*****     proton-induced knockout reaction calculation for exclusive processes     *****
!*****     triple- and quadruple-differential cross sections, analyzing powers,     *****
!*****         and momentum distributions in normal and inverse kinematics          *****
!*****                                                                              *****
!*****                                      by                                      *****
!*****                                                                              *****
!*****             Kazuyuki Ogata, Kazuki Yoshida, and Yoshiki Chazono              *****
!*****                                                                              *****
!****************************************************************************************

!=======================================================================
module consts
      implicit none
!=======================================================================
      real*8,parameter :: hc=197.3269788d0,ac=931.4940954d0,ec=137.035999139d0
      real*8,parameter :: chkcons=1.0d-2,chknorm=1.0d-3
      real*8,parameter :: pi=acos(-1.0d0),piad=pi/180.0d0
      real*8,parameter :: expmax=300.0d0,epslim=1.0d-8
      complex*16,parameter :: wi=(0.0d0,1.0d0)

end module consts
!=======================================================================



!=======================================================================
module dims
      implicit none
!=======================================================================
      integer,parameter :: nrdim=304
      integer,parameter :: nthdim=181,ldim=99,mdim=ldim
      integer,parameter :: ngdim=50
      integer,parameter :: nkapdim=51,nthnndim=181,nkfnndim=10
      integer,parameter :: nenndim= 401
      integer,parameter :: nepdim=82,nepfdim=201
      integer,parameter :: nzpldim=2001
      integer,parameter :: nq0dim=69,nthq0dim=13,nphq0dim=3
      integer,parameter :: nthqdim=13,nphqdim=3
      integer,parameter :: nkbdim=101,nkzdim=2*nkbdim-1
                           ! When nkbdim > 401, format #610&611 in "center"
                           ! must be edited accordingly.

end module dims
!=======================================================================



!=======================================================================
module array
      use dims
      implicit none
!=======================================================================
      integer lauto(3)

      real*8 ucou(3,nrdim)
      real*8 ylm0(nthdim,ldim,mdim)
      real*8 pltbl(ldim,nzpldim),dpltbl(ldim,nzpldim)
      real*8 dsigelm(nenndim,nthnndim)
      real*8 vpotw(nepdim,nrdim),wpotw(nepdim,nrdim)
      real*8 vlsw(nepdim,nrdim),wlsw(nepdim,nrdim)
      real*8 fnlwr(nepdim,nrdim),fnlwi(nepdim,nrdim)
      real*8 plcalsv(ngdim,ngdim,ldim)
      real*8 ddx(nkzdim,nkbdim)

      complex*16 fnl(3,nrdim),fnltbl(3,nepfdim,nrdim)
      complex*16 wf(3,ngdim,ngdim,ngdim),wfmtp(ngdim,ngdim,ngdim)
      complex*16 wfls(2,2,3,ngdim,ngdim,ngdim)
      complex*16 uopt(3,nrdim)
      complex*16 uls(3,nrdim)
      complex*16 uopttbl(3,nepfdim,nrdim),ulstbl(3,nepfdim,nrdim)
      complex*16 gma(nthqdim,nphqdim,nq0dim,nthq0dim,nphq0dim,16)

end module array
!=======================================================================



!=======================================================================
module angmom
      use dims
      implicit none
!=======================================================================
      real*8 faclog(ldim+mdim+1),dfaclog(mdim+1)
      integer memoam

end module angmom
!=======================================================================



!=======================================================================
module kibmod
      implicit none
!=======================================================================
      integer kibtbl,kibout,kibtmd,kiblg,kibpx,kibtr,kibtl,kibbs
      integer kibelm

end module kibmod
!=======================================================================



!=======================================================================
module nntbl
      use dims,only:nthq0dim,nphq0dim,nthqdim,nphqdim
      implicit none
!=======================================================================
      integer :: ifnntbl=1 ! flag for NN-table (gma) calculation in the 1st calc.
      integer ionsh
      integer nq0mn,nq0mx,nqmn,nqmx,nth0mn,nth0mx,nthmn,nthmx &
    &        ,nph0mn,nph0mx,nphmn,nphmx
      real*8 q0mn,q0mx,dq0,th0mn,qmn,qmx,dq,th0mx,thmn,thmx,dth &
    &       ,ph0mn,ph0mx,phmn,phmx,dph
      real*8 thq0mn,thq0mx,dthq0,phq0mn,phq0mx,dphq0 &
    &       ,thqmn,thqmx,dthq,phqmn,phqmx,dphq
      real*8 thq0tbl(nthq0dim),phq0tbl(nphq0dim)
      real*8 thqtbl(nthqdim),phqtbl(nphqdim)

end module nntbl
!=======================================================================



!=======================================================================
program pikoe1
      use kibmod,only:kibout
      implicit real*8(a-h,o-z)
!=======================================================================
      integer ipot(3)
      real*8 facv(3),facw(3),facvs(3),facws(3),beta(3)
      integer lmax(3),iedg(3),ims(3)
      integer icpot(3)

!-----------------------------------------------------------------------

!-------------------
!---  file open  ---
!-------------------

      call fopen


!-------------------------
!---  parameter input  ---
!-------------------------

      call input(limfs,ions,ifrm,imir,ical                   &
     &          ,zp,ap,za,aa,ikin,elab,ictrein               &
     &          ,ish,ebind,zsp,asp,betasp,ictrm              &
     &          ,fj,ls,sfac,nod,ibmc,rc,ictrc,a0c,rcl,ictrcl &
     &          ,ibms,v0ls,rs,ictrs,as,lmax                  &
     &          ,ivar,iex,fkncut,ixunt,kunt                  &
     &          ,ivvarl,varlmin,varlmax,dvarl                &
     &          ,ivthxl,thxlmin,thxlmax,dthxl                &
     &          ,ivphxl,phxlmin,phxlmax,dphxl                &
     &          ,ivet2l,et2lmin,et2lmax,det2l                &
     &          ,ivph2l,ph2lmin,ph2lmax,dph2l                &
     &          ,ielm,kinelm,ielmedg                         &
     &          ,rmax,dr,ngth,ngph,ngk1,ngph1,nrgmax         &
     &          ,ipot,facv,facw,facvs,facws,beta,ims,iedg    &
     &          ,detbl,itmdcal,icpot,ilscal)


!-------------------------
!---  tdx calculation  ---
!-------------------------

      call center(limfs,ions,ifrm,imir,ical,zp,ap,za,aa,betasp    &
     &           ,elab,ebind,zsp,asp,ibmc,rc,ictrc,a0c,rcl,ictrcl &
     &           ,ibms,v0ls,rs,ictrs,as,ielmedg,ims,iedg,ictrm    &
     &           ,iex,fkncut,ivar,ikin,ictrein                    &
     &           ,ivvarl,varlmin,varlmax,dvarl                    &
     &           ,ivthxl,thxlmin,thxlmax,dthxl                    &
     &           ,ivphxl,phxlmin,phxlmax,dphxl                    &
     &           ,ivet2l,et2lmin,et2lmax,det2l                    &
     &           ,ivph2l,ph2lmin,ph2lmax,dph2l                    &
     &           ,ielm,kinelm,fj,ls,sfac,nod,ish                  &
     &           ,rmax,dr,ngth,ngph,ngk1,ngph1,nrgmax             &
     &           ,ipot,facv,facw,facvs,facws,beta,lmax            &
     &           ,ixunt,kunt,detbl,icpot,itmdcal,ilscal)

      call clockm(it2)
      write(kibout,6020) it2*0.001
 6020 format(/1x, 45x,'>>> calculation completed (',f10.3,' sec)')

      stop 0

end program pikoe1
!=======================================================================



!=======================================================================
subroutine center(limfs,ions,ifrm,imir,ical,zp,ap,za,aa,betasp    &
     &           ,elab,ebind,zsp,asp,ibmc,rc,ictrc,a0c,rcl,ictrcl &
     &           ,ibms,v0ls,rs,ictrs,as,ielmedg,ims,iedg,ictrm    &
     &           ,iex,fkncut,ivar,ikin,ictrein                    &
     &           ,ivvarl,varlmin,varlmax,dvarl                    &
     &           ,ivthxl,thxlmin,thxlmax,dthxl                    &
     &           ,ivphxl,phxlmin,phxlmax,dphxl                    &
     &           ,ivet2l,et2lmin,et2lmax,det2l                    &
     &           ,ivph2l,ph2lmin,ph2lmax,dph2l                    &
     &           ,ielm,kinelm,fj,ls,sfac,nod,ish                  &
     &           ,rmax,dr,ngth,ngph,ngk1,ngph1,nrgmax             &
     &           ,ipot,facv,facw,facvs,facws,beta,lmax            &
     &           ,ixunt,kunt,detbl,icpot,itmdcal,ilscal)
      use consts
      use dims
      use array
      use kibmod
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 e1l(2),th1l(2),ph1l(2),e2l(2),th2l(2),fkbl(2),thbl(2),phbl(2)

      integer ipot(3)
      real*8 facv(3),facw(3),facvs(3),facws(3),beta(3)
      integer lmax(3),ims(3),iedg(3),ipread(3)
      integer iepfix(3),nepfmax(3)
      real*8 epfix(3),epfmin(3),r0cl(3)
      integer ictrcls(3)

      real*8 potv(nrdim),potw(nrdim),potd(nrdim),potrs(nrdim),potis(nrdim)

      real*8 rr(nrdim)

      real*8 vtot(nrdim),vcen(nrdim),vls(nrdim),vcou(nrdim)
      real*8 ffr(nrdim),ffrg(ngdim),fnlfac(nrdim),ffr_orig(nrdim)

      real*8 th(nthdim),costh(nthdim)

      real*8 ph1calg(ngdim),ph1weig(ngdim)
      real*8 fk1calg(ngdim),fk1weig(ngdim)
      real*8 thcalg(ngdim),thweig(ngdim)
      real*8 phcalg(ngdim),phweig(ngdim)
      real*8 rcalg(ngdim),rweig(ngdim)

      integer iapp(2)

      real*8 eelm(nenndim)
      real*8 rspsv(nrdim),wfsp(nrdim),fnlspw(nrdim),fnlsp(nrdim)

      real*8 epw(nepdim),rpw(nrdim)
      integer ir0cl(3)
      real*8 vcouw(nrdim)
      real*8 vpotw1d(nrdim),wpotw1d(nrdim)
      real*8 vlsw1d(nrdim),wlsw1d(nrdim)
      real*8 fnlwr1d(nrdim),fnlwi1d(nrdim)

      real*8 fkayi(3),thtkg(3),phikg(3),amasspi(3)
      real*8 pmi(3)
      integer icpot(3)

      real*8 ra(nrdim)
      character*12 sddxhd
      character*2 csunit

!-----------------------------------------------------------------------

!-----L-frame kinematics in the initial channel

      fm0=ap*ac
      fma=aa*ac
      fmsp=asp*ac

      if(ikin==0) then
         if(ictrein==0) then
           e0l=elab*nint(ap)+fm0
          else
           e0l=elab*ap+fm0
         end if
        fk0l=sqrt((e0l+fm0)*(e0l-fm0))/hc
        eal=fma
        fkal=0.0d0
       else
         if(ictrein==0) then
           eal=elab*nint(aa)+fma
          else
           eal=elab*aa+fma
         end if
        fkal=sqrt((eal+fma)*(eal-fma))/hc
        e0l=fm0
        fk0l=0.0d0
      end if

      betgl=hc*(fk0l+fkal)/(eal+e0l)
      gamgl=1.0d0/sqrt(1.0d0-betgl**2)
      betlg=-betgl
      gamlg=gamgl

!-----V-frame kinematics in the initial channel

      if(ikin==0) then
        betvl=hc*fk0l/e0l
        gamvl=1.0d0/sqrt(1.0d0-betvl**2)
        e0v=fm0
        call ltrz(1,eal,1,0.0d0,0.0d0,betvl,gamvl,fma,eav,fkav,thav,phav)
        fk0v=0.0d0
        fk0vz=0.0d0
        fkavz=-fkav
        betgv=-hc*fkav/(e0v+eav)
       else
        betvl=hc*fkal/eal
        gamvl=1.0d0/sqrt(1.0d0-betvl**2)
        call ltrz(1,e0l,1,0.0d0,0.0d0,betvl,gamvl,fm0,e0v,fk0v,th0v,ph0v)
        eav=fma
        fk0vz=-fk0v
        fkav=0.0d0
        fkavz=0.0d0
        betgv=-hc*fk0v/(e0v+eav)
      end if
      gamgv=1.0d0/sqrt(1.0d0-betgv**2)
      betvg=-betgv
      gamvg=gamgv

!-----A-frame kinematics in the initial channel

      if(ikin==0) then
        fk0a=fk0l
        fk0az=fk0a
        e0a=e0l
        betal=0.0d0
        betga=betgl
       else
        e0a=e0v
        eaa=eav
        fk0a=fk0v
        fk0az=-fk0a
        betal=betvl
        betga=betgv
      end if
      t0a=e0a-fm0
      gamal=1.0d0/sqrt(1.0d0-betal**2)
      gamga=1.0d0/sqrt(1.0d0-betga**2)
      betag=-betga
      gamag=gamga


!-------------------------------------------------
!---  bound state wave function (radial part)  ---
!-------------------------------------------------

      nrmax=nint(rmax/dr)+1
      if(nrmax>nrdim) then
        write(*,*) 'ERROR: nrmax > nrdim'
        write(*,'(a,i0,a,i0)') ' nrmax=',nrmax,', nrdim=',nrdim
        stop 1
      end if
      do ir=1,nrmax
       rr(ir)=(ir-1)*dr
      end do

      izsp=nint(zsp)
      zb=za-zsp
      zzeb=zsp*zb*hc/ec

      call potbs(izsp,asp,za,aa,fj,ls,ebind,zzeb,nrmax,dr           &
     &          ,ibmc,rc,ictrc,a0c,rcl,ictrcl,ibms,v0ls,rs,ictrs,as &
     &          ,radi,rads,radc,vtot,vcen,vls,vcou,vdepc,vdeps)

      pmass=asp
      if(ictrm==1) then
        ffmm=2.0d0*asp*ac/hc**2
       else if(ictrm==2) then
        fmb=fma-fmsp  ! note: fmb actually depends on ebind.
        fmusp=asp*ac*fmb/(asp*ac+fmb)
        ffmm=2.0d0*fmusp/hc**2
      end if

      if(ish<10) then
        amass=aa-asp
        pvc=vdepc
        if(ish==1) pvc=-ebind
        pvs=vdeps
        ish0=ish

        call bound2(ffmm,zzeb,ls,nrmax,dr,energy,ffr,vtot, &
     &              vcen,vls,vcou,radi,radc,pmass,         &
     &              pvc,pvs,ish,nod,vdeptho,wlso)

        if(ish0==0) ebind=-energy

       else

        read(ish,*) rspmax,drsp
        if(rspmax<=0.0d0) then
          write(*,*) 'ERROR: rspmax must be positive'
          write(*,'(a,f0.4)') ' rspmax=',rspmax
          stop 1
        end if

        if(drsp<=0.0d0) then
          write(*,*) 'ERROR: drsp must be positive'
          write(*,'(a,f0.4)') ' drsp=',drsp
          stop 1
        end if
        nrspmax=nint(rspmax/drsp)+1
        iflgnlsp=0
        do irsp=1,nrspmax
         read(ish,502) rspsv(irsp),wfsp(irsp),fnlspw(irsp)
  502    format(f10.0,2e20.12)
         if(abs(fnlspw(irsp))>epslim) iflgnlsp=1
        end do
        if(iflgnlsp==0) fnlspw(:)=1.0d0
        do ir=1,nrmax
         r=(ir-1)*dr
         ffr(ir)=suphod(r,rspsv,drsp,wfsp,nrspmax,nrdim)
         fnlsp(ir)=suphod(r,rspsv,drsp,fnlspw,nrspmax,nrdim)
        end do
      end if


!--------------------------------------------------------------
!---  nonlocal correction to the bound-state wave function  ---
!--------------------------------------------------------------

      fnlfac(:)=1.0d0
      ffr_orig(:)=ffr(:)
      if(abs(betasp)>epslim) then
        facbtsp=ffmm*betasp**2/4.0d0
        sumsp=0.0d0
        do ir=1,nrmax
         if(betasp>epslim) then
           irc=ir
!-----     vls(1) may have a "core" that causes a numerical problem.
           if(ir==1) irc=2
           fnlfac(ir)=(1.0d0-facbtsp*(vtot(irc)-vcou(irc)))**(-0.5d0)
          else
           fnlfac(ir)=fnlsp(ir)
         end if
         ffr(ir)=ffr(ir)*fnlfac(ir)
         sumsp=sumsp+ffr(ir)**2*rr(ir)**2
        end do
        sumsp=sumsp*dr
        do ir=1,nrmax
         ffr(ir)=ffr(ir)/sqrt(sumsp)
        end do
      end if

      if(kibbs>0) then 
        write(kibbs,620)
  620   format(4x,'r[fm]',10x,'wf[fm-3/2]',5x,'NLC',12x,'corrected-wf[fm-3/2]')
        do ir=1,nrmax
         r=(ir-1)*dr
         write(kibbs,'(4e15.5)') r,ffr_orig(ir),fnlfac(ir),ffr(ir)
        end do
      end if

      do ir=1,nrmax
       ra(ir)=(ir-1)*dr
      end do
      call glweight(nrgmax,0.d0,rmax,rcalg,rweig)
      do ir=1,nrgmax
       r=rcalg(ir)
       ffrg(ir)=polintm(ra(1:nrmax),ffr(1:nrmax),nrmax,r,3)
      end do

!-----particle masses in the final channel

      fmsp=asp*ac
      fmb=fma-fmsp+ebind
      ab=fmb/ac
      if(iex==1) then
        fm1=fmsp
        z1=zsp
        fm2=fm0
        z2=zp
       else
        fm1=fm0
        z1=zp
        fm2=fmsp
        z2=zsp
      end if

!-----G-frame kinematics of particle 0 and a

      call ltrz(2,fk0l,1,0.0d0,0.0d0,betgl,gamgl,fm0,e0,fk0,th0,ph0)
      call vecxyz(1,fk0,th0,ph0, fk0x,fk0y,fk0z)
      call ltrz(2,fkal,1,0.0d0,0.0d0,betgl,gamgl,fma,ea,fka,tha,pha)
      call vecxyz(1,fka,tha,pha, fkax,fkay,fkaz)

!-----output of the kinematics of particles 0 and a

      write(kibout,601) e0l,e0l-fm0,fk0l, eal,eal-fma,fkal     &
     &                 ,e0,e0-fm0,fk0,fk0z, ea,ea-fma,fka,fkaz &
     &                 ,e0v,e0v-fm0,fk0v, eav,eav-fma,fkav
  601 format(/3x,'-- initial-channel kinematics outputs --'     &
     &       /8x,'e0l t0l  fk0l                :',2f10.2,f10.5  &
     &       /8x,'eal tal  fkal                :',2f10.2,f10.5  &
     &       /8x,'e0  t0   fk0  fk0z           :',2f10.2,2f10.5 &
     &       /8x,'ea  ta   fka  fkaz           :',2f10.2,2f10.5 &
     &       /8x,'e0v t0p  fk0v                :',2f10.2,f10.5  &
     &       /8x,'eav tap  fkav                :',2f10.2,f10.5)

      write(kibout,602) betgl,gamgl,betvl,gamvl,betvg,gamvg
  602 format(/3x,'-- Lorentz factors  --'                 &
     &       /8x,'betgl  gamgl  (g <-- l)      :',2e20.12 &
     &       /8x,'betvl  gamvl  (v <-- l)      :',2e20.12 &
     &       /8x,'betvg  gamvg  (v <-- g)      :',2e20.12)


!--------------------------------------
!---  distorting potential read-in  ---
!--------------------------------------

      ir0cl(1:3)=0
      ictrcls(1:3)=1
      ipread(1:3)=0
      nepfmax(1:3)=0

      if(ipot(1)>9) then
        ipread(1)=1
        call potread(ipot(1),1,nrmax,dr                            &
     &              ,vpotw1d,wpotw1d,vlsw1d,wlsw1d,epw,rpw         &
     &              ,vpotw,wpotw,vcouw,vlsw,wlsw,ir0cl,r0cl(1)     &
     &              ,ictrcls(1)                                    &
     &              ,iepfix(1),epfix(1),epfmin(1),nepfmax(1),detbl &
     &              ,uopt,uls,ucou,uopttbl,ulstbl                  &
     &              ,fnlwr,fnlwi,fnlwr1d,fnlwi1d,fnl,fnltbl)
      end if

      if(ipot(2)>9) then
        if(ipot(2)==ipot(1)) then
          iepfix(2) =iepfix(1)
          epfix(2)  =epfix(1)
          epfmin(2) =epfmin(1)
          nepfmax(2)=nepfmax(1)
          ir0cl(2)  =ir0cl(1)
          r0cl(2)   =r0cl(1)
          ictrcls(2)=ictrcls(1)
          if(ir0cl(2)==2) then
            if(abs(z1*zb-zp*za)>epslim) then
              write(*,'(a)') 'ERROR: charge product changed but Coulomb potential not'
              stop 1
            end if
            do ir=1,nrmax
             ucou(2,ir)=ucou(1,ir)
            end do
          end if
          if(iepfix(2)==1) then
            do ir=1,nrmax
             uopt(2,ir)=dble(uopt(1,ir))*facv(2) +imag(uopt(1,ir))*facw(2)*wi
             uls (2,ir)=dble(uls (1,ir))*facvs(2)+imag(uls(1,ir))*facws(2)*wi
             fnl(2,ir)=fnl(1,ir)
            end do
           else
            do iepf=1,nepfmax(2)
             do ir=1,nrmax
              uopttbl(2,iepf,ir)=dble(uopttbl(1,iepf,ir))*facv(2) &
     &                          +imag(uopttbl(1,iepf,ir))*facw(2)*wi
              ulstbl (2,iepf,ir)=dble(ulstbl(1,iepf,ir))*facvs(2) &
     &                          +imag(ulstbl(1,iepf,ir))*facws(2)*wi
              fnltbl(2,iepf,ir)=fnltbl(1,iepf,ir)
             end do
            end do
          end if
         else
          ipread(2)=1
          call potread(ipot(2),2,nrmax,dr                            &
     &                ,vpotw1d,wpotw1d,vlsw1d,wlsw1d,epw,rpw         &
     &                ,vpotw,wpotw,vcouw,vlsw,wlsw,ir0cl,r0cl(2)     &
     &                ,ictrcls(2)                                    &
     &                ,iepfix(2),epfix(2),epfmin(2),nepfmax(2),detbl &
     &                ,uopt,uls,ucou,uopttbl,ulstbl                  &
     &                ,fnlwr,fnlwi,fnlwr1d,fnlwi1d,fnl,fnltbl)
        end if
      end if

      if(ipot(3)>9) then
        if(ipot(3)==ipot(1)) then
          iepfix(3) =iepfix(1)
          epfix(3)  =epfix(1)
          epfmin(3) =epfmin(1)
          nepfmax(3)=nepfmax(1)
          ir0cl(3)  =ir0cl(1)
          r0cl(3)   =r0cl(1)
          ictrcls(3)=ictrcls(1)
          if(ir0cl(3)==2) then
            if(abs(z2*zb-zp*za)>epslim) then
              write(*,'(a)') 'ERROR: charge product changed but Coulomb potential not'
              stop 1
            end if
            do ir=1,nrmax
             ucou(3,ir)=ucou(1,ir)
            end do
          end if
          if(iepfix(3)==1) then
            do ir=1,nrmax
             uopt(3,ir)=dble(uopt(1,ir))*facv(3)+imag(uopt(1,ir))*facw(3)*wi
             uls (3,ir)=dble(uls (1,ir))*facvs(3)+imag(uls (1,ir))*facws(3)*wi
             fnl(3,ir)=fnl(1,ir)
            end do
           else
            do iepf=1,nepfmax(3)
             do ir=1,nrmax
              uopttbl(3,iepf,ir)=dble(uopttbl(1,iepf,ir))*facv(3) &
     &                          +imag(uopttbl(1,iepf,ir))*facw(3)*wi
              ulstbl(3,iepf,ir)= dble(ulstbl(1,iepf,ir))*facvs(3) &
     &                          +imag(ulstbl(1,iepf,ir))*facws(3)*wi
              fnltbl(3,iepf,ir)=fnltbl(1,iepf,ir)
             end do
            end do
          end if
         else if(ipot(3)==ipot(2)) then
          iepfix(3) =iepfix(2)
          epfix(3)  =epfix(2)
          epfmin(3) =epfmin(2)
          nepfmax(3)=nepfmax(2)
          ir0cl(3)  =ir0cl(2)
          r0cl(3)   =r0cl(2)
          ictrcls(3)=ictrcls(2)
          if(ir0cl(3)==2) then
            if(abs(z2*zb-z1*zb)>epslim) then
              write(*,'(a)') 'ERROR: charge product changed but Coulomb potential not'
              stop 1
            end if
            do ir=1,nrmax
             ucou(3,ir)=ucou(2,ir)
            end do
          end if
          if(iepfix(3)==1) then
            do ir=1,nrmax
             uopt(3,ir)=dble(uopt(2,ir))*facv(3) +imag(uopt(2,ir))*facw(3)*wi
             uls (3,ir)=dble(uls (2,ir))*facvs(3)+imag(uls(2,ir))*facws(3)*wi
             fnl(3,ir)=fnl(2,ir)
            end do
           else
            do iepf=1,nepfmax(3)
             do ir=1,nrmax
              uopttbl(3,iepf,ir)=dble(uopttbl(2,iepf,ir))*facv(3) &
     &                          +imag(uopttbl(2,iepf,ir))*facw(3)*wi
              ulstbl (3,iepf,ir)=dble(ulstbl(2,iepf,ir))*facvs(3) &
     &                          +imag(ulstbl(2,iepf,ir))*facws(3)*wi
              fnltbl(3,iepf,ir)=fnltbl(2,iepf,ir)
             end do
            end do
          end if
         else
          ipread(3)=1
          call potread(ipot(3),3,nrmax,dr                            &
     &                ,vpotw1d,wpotw1d,vlsw1d,wlsw1d,epw,rpw         &
     &                ,vpotw,wpotw,vcouw,vlsw,wlsw,ir0cl,r0cl(3)     &
     &                ,ictrcls(3)                                    &
     &                ,iepfix(3),epfix(3),epfmin(3),nepfmax(3),detbl &
     &                ,uopt,uls,ucou,uopttbl,ulstbl                  &
     &                ,fnlwr,fnlwi,fnlwr1d,fnlwi1d,fnl,fnltbl)
        end if
      end if

      do iptcl=1,3
       if(nepfmax(iptcl)>nepfdim) then
         write(*,*) 'ERROR: nepfmax > nepfdim'
         write(*,'(a,i0,a,i0)') ' nepfmax=',nepfmax(iptcl),', nepfdim=',nepfdim
         stop 1
       end if
       if(ipread(iptcl)==0) cycle
       if(iepfix(iptcl)==1) then
         do ir=1,nrmax
          uopt(iptcl,ir)=dble(uopt(iptcl,ir))*facv(iptcl) &
     &                  +imag(uopt(iptcl,ir))*facw(iptcl)*wi
          uls (iptcl,ir)=dble(uls(iptcl,ir))*facvs(iptcl) &
     &                  +imag(uls(iptcl,ir))*facws(iptcl)*wi
         end do
        else
         do iepf=1,nepfmax(iptcl)
          do ir=1,nrmax
           uopttbl(iptcl,iepf,ir)=dble(uopttbl(iptcl,iepf,ir))*facv(iptcl) &
     &                           +imag(uopttbl(iptcl,iepf,ir))*facw(iptcl)*wi
           ulstbl (iptcl,iepf,ir)=dble(ulstbl(iptcl,iepf,ir))*facvs(iptcl) &
     &                           +imag(ulstbl(iptcl,iepf,ir))*facws(iptcl)*wi
          end do
         end do
       end if
      end do


!---------------------------
!---  coulomb potential  ---
!---------------------------

      do iptcl=1,3
       afac=ab
       if(iptcl==1) afac=aa
       if(iptcl==1) then
         amsct=ap
         amtgt=aa
        else if(iptcl==2) then
         amsct=fm1/ac
         amtgt=ab
        else
         amsct=fm2/ac
         amtgt=ab
       end if
       if(ictrcls(iptcl)==1) then
         faccl=amtgt**(1.0d0/3.0d0)
        else if(ictrcls(iptcl)==2) then
         faccl=amsct**(1.0d0/3.0d0)+amtgt**(1.0d0/3.0d0)
        else if(ictrcls(iptcl)==3) then
         faccl=1.0d0
        else
         write(*,*) 'ERROR: ictrcls must be 1--3'
         write(*,'(a,i0,a,i0)') ' ictrcls=',ictrcls(iptcl),', for particle=',iptcl-1
         stop 1
       end if
       if(ir0cl(iptcl)==0) then
!-----   kd parameter
         rc=(1.198d0+0.697d0*afac**(-2.0d0/3.0d0) &
     &        +12.994d0*afac**(-5.0d0/3.0d0))*afac**(1.0d0/3.0d0)
        else if(ir0cl(iptcl)==1) then
         rc=r0cl(iptcl)*faccl
        else if(ir0cl(iptcl)==2) then
         cycle
       end if
       if(iptcl==1) then
         zze=zp*za*hc/ec
       else if(iptcl==2) then
         zze=z1*zb*hc/ec
       else if(iptcl==3) then
         zze=z2*zb*hc/ec
       end if
       do ir=1,nrmax
        rcal=rr(ir)
        if(ir==1) rcal=rcal+1.0d-10
        ucou(iptcl,ir)=zze/rcal
        if(rcal<=rc) ucou(iptcl,ir)=zze/(2.0d0*rc)*(3.0d0-rcal**2/rc**2)
       end do
      end do


!-----------------------------------------------
!---  three-body kinematics and observables  ---
!-----------------------------------------------

      fslim=dble(limfs)

      nvarlmax=1
      nthxlmax=1
      nphxlmax=1
      net2lmax=1
      nph2lmax=1

      i1st=1
      io=0

      call header(ivar,ifrm,imir,ical,ixunt,kunt)

      csunit='mb'
      if(ixunt==1) csunit='ub'

      if(kunt==0) then
        fkfac=1.0d0
       else if(kunt==1) then
        fkfac=hc
       else if(kunt==2) then
        fkfac=hc*1.0d-3
      end if

      if(ivar<9) then
        if(ivvarl/=0) nvarlmax=nint((varlmax-varlmin)/dvarl)+1
        if(ivthxl/=0) nthxlmax=nint((thxlmax-thxlmin)/dthxl)+1
        if(ivphxl/=0) nphxlmax=nint((phxlmax-phxlmin)/dphxl)+1
        if(ivet2l/=0) net2lmax=nint((et2lmax-et2lmin)/det2l)+1
        if(ivph2l/=0) nph2lmax=nint((ph2lmax-ph2lmin)/dph2l)+1

        dthxlrad=dthxl*piad
        dphxlrad=dphxl*piad
        if(ivar==3) then
          de2l=det2l
         else
          dth2lrad=det2l*piad
        end if
        dph2lrad=dph2l*piad

        if(ivvarl==0) dvarl   =1.d0
        if(ivthxl==0) dthxlrad=1.d0
        if(ivphxl==0) dphxlrad=1.d0
        if(ivet2l==0) then
          if(ivar==3) then
            de2l=1.0d0
           else
            dth2lrad=1.0d0
          end if
        end if
        if(ivph2l==0) dph2lrad=1.d0

        facid=1.0d0
        if(nint(ap)==nint(asp) .and. nint(zp)==nint(zsp)) facid=0.5d0
        tottdx=0.d0

        do iv=1,nvarlmax
         varl=(iv-1)*dvarl+varlmin

         do ithx=1,nthxlmax
          thxl=(ithx-1)*dthxl+thxlmin
          do iphx=1,nphxlmax
           phxl=(iphx-1)*dphxl+phxlmin
           do iet2=1,net2lmax
            et2l=(iet2-1)*det2l+et2lmin
            do iph2=1,nph2lmax
             ph2l=(iph2-1)*dph2l+ph2lmin

             tdx=0.0d0
             ay=0.0d0

             call kinema3bl(e0l,fk0l,eal,fkal,fmb,fm1,fm2 &
     &                     ,ivar,varl,thxl,phxl,et2l,ph2l &
     &                     ,e1l,th1l,ph1l,e2l,th2l,fkbl,thbl,phbl,nsol,iapp)

             if(nsol==0) then
               if(ions==0) cycle
               call tblout(e1l(1),th1l(1),ph1l(1),e2l(1),th2l(1),ph2l       &
     &                    ,fkbl(1),thbl(1),phbl(1)                          &
     &                    ,ikin,ifrm,imir,ivar,nsol,0,fm1,fm2,fmb,kunt,ions &
     &                    ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
               io=io+1
               if(io*1.22d-4>fslim) then
                 write(kibout,690)
  690            format(/1x,'  *** file size exeeds fslim (mb): program terminated ***')
                 go to 1000
               end if
               cycle
             end if

             do isol=1,2
              if(iapp(isol)==0) cycle

              io=io+1

              call ltrz(1,e1l(isol), 1,th1l(isol),ph1l(isol),betgl,gamgl &
     &                 ,fm1,  e1,fk1,th1,ph1)
              call ltrz(1,e2l(isol), 1,th2l(isol),ph2l,betgl,gamgl       &
     &                 ,fm2,  e2,fk2,th2,ph2)
              call ltrz(2,fkbl(isol),1,thbl(isol),phbl(isol),betgl,gamgl &
     &                 ,fmb,  eb,fkb,thb,phb)

              fk1l=sqrt((e1l(isol)+fm1)*(e1l(isol)-fm1))/hc
              fk2l=sqrt((e2l(isol)+fm2)*(e2l(isol)-fm2))/hc

              if(ikin==0) then
                t1a=e1l(isol)-fm1
                t2a=e2l(isol)-fm2
               else
                call ltrz(1,e1, 1,th1,ph1,betag,gamag,fm1,e1a,fk1a,th1a,ph1a)
                call ltrz(1,e2, 1,th2,ph2,betag,gamag,fm2,e2a,fk2a,th2a,ph2a)
                t1a=e1a-fm1
                t2a=e2a-fm2
              end if

              call vecxyz(1,fkb,thb,phb, fkbx,fkby,fkbz)
              fknz=-fkbz-(aa-asp)/aa*fk0z
              fknb=sqrt(fkbx**2+fkby**2)
              fkn=sqrt(fknz**2+fknb**2)
              if(fkncut>epslim .and. fkn>fkncut) then
                if(ions==0) then
                  io=io-1
                  cycle
                end if
                call tblout(e1l(isol),th1l(isol),ph1l(isol)           &
     &                     ,e2l(isol),th2l(isol),ph2l                 &
     &                     ,fkbl(isol),thbl(isol),phbl(isol)          &
     &                     ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb &
     &                     ,kunt,ions                                 &
     &                     ,betgl,gamgl,betvl,gamvl,1,ical,tdx,ay)
                if(io*1.22d-4>fslim) then
                  write(kibout,690)
                  go to 1000
                end if
                cycle
              end if

              if(ical==0) then
                call tblout(e1l(isol),th1l(isol),ph1l(isol)           &
     &                     ,e2l(isol),th2l(isol),ph2l                 &
     &                     ,fkbl(isol),thbl(isol),phbl(isol)          &
     &                     ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb &
     &                     ,kunt,ions                                 &
     &                     ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
                if(io*1.22d-4>fslim) then
                  write(kibout,690)
                  go to 1000
                end if
                cycle
              end if

              call tdxqm(rmax,dr,fkayi,thtkg,phikg                    &
     &                  ,elab,amasspi,e0,e1,e2,ea,eb,t0a,t1a,t2a,ylm0 &
     &                  ,zp,z1,z2,za,zb,pwsum3                        &
     &                  ,rcn,ircn,pmi,beta,sfac,ucou,dthnn,dsigelm    &
     &                  ,wf,wfmtp,ffrg,nrmax,i1st,fj,ls,lmax,icpot    &
     &                  ,ngth,ngph,ielm,kinelm,sig                    &
     &                  ,ipot,facv,facw,facvs,facws,eelm,iex,io       &
     &                  ,izsp,ap,asp,aa,ab,th,costh                   &
     &                  ,fnl,fnltbl,ims,iedg,ielmedg                  &
     &                  ,iepfix,epfix,nepfmax,epfmin                  &
     &                  ,detbl,uopttbl,ulstbl,uopt                    &
     &                  ,potv,potw,potd,potrs,potis                   &
     &                  ,fm0,fm1,fm2,fmsp                             &
     &                  ,nthnnmax,thnnmin,nennmax,sigelm              &
     &                  ,fk0,th0,ph0,fk1,th1,ph1,fk2,th2,ph2          &
     &                  ,thcalg,thweig,nthgmax                        &
     &                  ,phcalg,phweig,nphgmax                        &
     &                  ,rcalg,rweig,nrgmax                           &
     &                  ,pltbl,dpltbl,dzpl,nzplmax                    &
     &                  ,hv0,c0r,fmuelm,const,tdxr                    &
     &                  ,itmdcal,ilscal,wfls,uls,ay)

              phvol=0.0d0
              fjacob=1.0d0
              if(ifrm==0) then
                fktotiz=fk0l+fkal
                ebl=sqrt((hc*fkbl(isol))**2+fmb**2)
                phvol=fkin(fk1l,e1l(isol),th1l(isol),ph1l(isol),fk2l,e2l(isol)  &
     &                    ,th2l(isol),ph2l,fkbl(isol),ebl,thbl(isol),phbl(isol) &
     &                    ,fktotiz,ivar)
                fjacob=fjac(e1,e2,eb,e1l(isol),e2l(isol),ebl)
               else if(ifrm==1) then
                phvol=fkin(fk1,e1,th1,ph1,fk2,e2,th2,ph2,fkb,eb,thb,phb,0.0d0,ivar)
               else if(ifrm==2) then
                call ltrz(1,e1l(isol), 1,th1l(isol),ph1l(isol),betvl,gamvl &
     &                   ,fm1,  e1v,fk1v,th1v,ph1v)
                call ltrz(1,e2l(isol), 1,th2l(isol),ph2l,betvl,gamvl       &
     &                   ,fm2,  e2v,fk2v,th2v,ph2v)
                call ltrz(2,fkbl(isol),1,thbl(isol),phbl(isol),betvl,gamvl &
     &                   ,fmb,  ebv,fkbv,thbv,phbv)
                fktotiz=-fk0v-fkav
                phvol=fkin(fk1v,e1v,th1v,ph1v,fk2v,e2v,th2v,ph2v &
     &                    ,fkbv,ebv,thbv,phbv,fktotiz,ivar)
                fjacob=fjac(e1,e2,eb,e1v,e2v,ebv)
              end if

              tdx=phvol*fjacob*tdxr
              if(ixunt==1) tdx=tdx*1.0d3
              if(ifrm==0) then
                tdxinteg=tdx
                if(nthxlmax>1) tdxinteg=tdxinteg*sin(th1l(isol)*piad)
                if(net2lmax>1 .and. ivar/=3) tdxinteg=tdxinteg*sin(th2l(isol)*piad)
                if(nvarlmax>1.and.ivar==2) tdxinteg=tdxinteg*varl**2
                if(nvarlmax>1.and.(iv == 1.or.iv == nvarlmax)) tdxinteg=tdxinteg/2.d0
                if(ivthxl==1 .and.(ithx==1.or.ithx==nthxlmax)) tdxinteg=tdxinteg/2.d0
                if(ivphxl==1 .and.(iphx==1.or.iphx==nphxlmax)) tdxinteg=tdxinteg/2.d0
                if(ivet2l==1 .and.(iet2==1.or.iet2==net2lmax)) tdxinteg=tdxinteg/2.d0
                if(ivph2l==1 .and.(iph2==1.or.iph2==nph2lmax)) tdxinteg=tdxinteg/2.d0
                tottdx=tottdx+tdxinteg
              end if

              call tblout(e1l(isol),th1l(isol),ph1l(isol)           &
     &                   ,e2l(isol),th2l(isol),ph2l                 &
     &                   ,fkbl(isol),thbl(isol),phbl(isol)          &
     &                   ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb &
     &                   ,kunt,ions                                 &
     &                   ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
              if(io*1.22d-4>fslim) then
                write(kibout,690)
                go to 1000
              end if

             end do ! isol

            end do ! iph2
           end do ! iet2
          end do ! iphx
         end do ! ithx
        end do ! iv

        facint=1.0d0
        if(ivphxl==0 .and. ivvarl*ivthxl*ivet2l*ivph2l/=0) facint=2.0d0*pi
        if(ifrm==0) then
          if(ivar==3) then
            tottdx=tottdx*dvarl*dthxlrad*dphxlrad*de2l*dph2lrad*facid
           else
            tottdx=tottdx*dvarl*dthxlrad*dphxlrad*dth2lrad*dph2lrad*facid
          end if
          if(nint(facint)==1) then
            write(kibout,653) tottdx,csunit
           else
            write(kibout,654) tottdx,csunit,tottdx*facint,csunit
          end if
  653     format(/3x,'-- integrated value of the calculated TDX --' &
     &           /8x,f12.2,1x,a2)
  654     format(/3x,'-- integrated value of the calculated TDX --' &
     &           /8x,f12.2,1x,a2,3x,'(x 2pi =',f12.2,1x,a2,')')
          if(nint(ap)==1 .and. nint(asp)==1) then
            call signn(elab,sigpp,sigpn)
            signntot=sigpp
            if(abs(zp-zsp)>epslim) signntot=sigpn
            write(kibout,655) signntot
  655       format(6x,'(cf. NN total cross section:',f7.2,' mb)')
          end if
        end if

       else if(ivar==9) then

        nph1gmax=ngph1
        call glweight(nph1gmax,0.d0,2.d0*pi,ph1calg,ph1weig)

        if(ivvarl/=0) nvarlmax=nint((varlmax-varlmin)/dvarl)+1
        if(ivthxl/=0) nthxlmax=nint(thxlmax/dthxl)+1

        fkbazmin=varlmin
        fkbazmax=varlmax
        dkbaz=dvarl
        fkbabmin=thxlmin
        fkbabmax=thxlmax
        dkbab=dthxl
        nkbazmax=nvarlmax
        nkbabmax=nthxlmax
        if(nkbazmax>nkzdim) then
          write(*,*) 'ERROR: nkbazmax > nkzdim'
          write(*,'(a,i0,a,i0)') ' nkbazmax=',nkbazmax,', nkzdim=',nkzdim
          stop 1
        end if
        if(nkbabmax>nkbdim) then
          write(*,*) 'ERROR: nkbabmax > nkbdim'
          write(*,'(a,i0,a,i0)') ' nkbabmax=',nkbabmax,', nkbdim=',nkbdim
          stop 1
        end if

        facid=1.0d0
        if(nint(ap)==nint(asp) .and. nint(zp)==nint(zsp)) facid=0.5d0
        sumsig=0.0d0
        ddx(:,:)=0.0d0
        if(kunt==0) then
          if(imir==0) then
            sddxhd=' kbaz-kbab  '
           else
            sddxhd=' kbaz:m-kbab'
          end if
         else
          if(imir==0) then
            sddxhd=' pbaz-pbab  '
           else
            sddxhd=' pbaz:m-pbab'
          end if
        end if
        write(kibtbl,610) sddxhd,(((i-1)*dkbab+fkbabmin)*fkfac,i=1,nkbabmax)
  610   format(a12,401f13.5)

        do ikbaz=1,nkbazmax
         fkbaz=(ikbaz-1)*dkbaz+fkbazmin
         fkbazo=fkbaz
         if(imir==1) fkbazo=-fkbazo
         sumkbb=0.0d0
         do ikbab=1,nkbabmax
          fkbab=(ikbab-1)*dkbab
          if(fkbab<fkbabmin-epslim) cycle
          fkba=sqrt(fkbab**2+fkbaz**2)
          eba=sqrt(hc**2*fkba**2+fmb**2)

          ay=0.0d0

          fkbz=gamga*(fkbaz-betga*eba/hc)
          fkbx=fkbab
          fkby=0.0d0
          fkb=sqrt(fkbx**2+fkby**2+fkbz**2)
          eb=gamga*(eba-betga*fkbaz*hc)

          fmolb=eb/eba

          fknz=-fkbz-(aa-asp)/aa*fk0z
          fknx=-fkbab
          fknb=fkbab
          fkn=sqrt(fknz**2+fknx**2)
          if(fkncut>epslim .and. fkn>fkncut) cycle

          qaz=fk0az-fkbaz
          qax=-fkbab
          qa2=qaz**2+qax**2
          qa=sqrt(qa2)
          cosqa=qaz/qa
          sinqa=sqrt(1.0d0-cosqa**2)
          epsbara=(e0a+fma-eba)/hc

          faca=fm1/hc
          facb=(qa2-epsbara**2+(fm2/hc)**2-(fm1/hc)**2)/(2.0d0*epsbara)
          facc=epsbara/qa
          sqfac2=facb**2*facc**2+faca**2-faca**2*facc**2

          fm12heavy=max(fm1,fm2)
          e1e2summin=sqrt((qa*hc)**2+fm12heavy**2)
          if(epsbara*hc<e1e2summin) cycle
          if(facb>0.0d0) then
            write(*,*) 'ERROR: facb must be negative'
            write(*,'(a,e13.5)') 'facb=',facb
            stop 1
          end if
          if(facb<=-faca) then
            sqfac=sqrt(sqfac2)
            fk1amin=( facb*facc+facc*sqfac)/(facc**2-1.0d0)
            fk1amax=(-facb*facc+facc*sqfac)/(facc**2-1.0d0)
           else
            if(sqfac2<0.020) then
              sumk1a=0.0d0
              cycle
            end if
            sqfac=sqrt(sqfac2)
            fk1amin=(-facb*facc-facc*sqfac)/(facc**2-1.0d0)
            fk1amax=(-facb*facc+facc*sqfac)/(facc**2-1.0d0)
          end if

          nk1gmax=ngk1
          call glweight(nk1gmax,fk1amin,fk1amax,fk1calg,fk1weig)

          sumk1a=0.0d0
          do ik1ag=1,nk1gmax
           fk1a=fk1calg(ik1ag)
           e1a=sqrt(hc**2*fk1a**2+fm1**2)
           t1a=e1a-fm1
           e2a=epsbara*hc-e1a
           fk2a=sqrt((e2a+fm2)*(e2a-fm2))/hc
           t2a=e2a-fm2
           if(fk1a<1.0d-10) cycle
           cos1qa=facc*(sqrt(fk1a**2+faca**2)+facb)/fk1a
           sin1qa=sqrt(1.0d0-cos1qa**2)
           sumph1a=0.0d0

           do iph1g=1,nph1gmax
            ph1a=ph1calg(iph1g)
            if(nph1gmax==1) ph1a=0.0d0
            fk1ax=fk1a*(sin1qa*cos(ph1a)*cosqa-cos1qa*sinqa)
            fk1ay=fk1a*sin1qa*sin(ph1a)
            fk1az=fk1a*(sin1qa*cos(ph1a)*sinqa+cos1qa*cosqa)

            fk1x=fk1ax
            fk1y=fk1ay
            fk1z=gamga*(fk1az-betga*e1a/hc)
            call vecthph(1,fk1x,fk1y,fk1z ,fk1,th1,ph1)
            e1=sqrt(hc**2*fk1**2+fm1**2)

            fk2x=-fk1x-fkbx
            fk2y=-fk1y
            fk2z=-fk1z-fkbz
            call vecthph(1,fk2x,fk2y,fk2z ,fk2,th2,ph2)
            e2=sqrt(hc**2*fk2**2+fm2**2)

            fmol=fmolb*e1*e2/(e1a*e2a)

            call tdxqm(rmax,dr,fkayi,thtkg,phikg                    &
     &                ,elab,amasspi,e0,e1,e2,ea,eb,t0a,t1a,t2a,ylm0 &
     &                ,zp,z1,z2,za,zb,pwsum3                        &
     &                ,rcn,ircn,pmi,beta,sfac,ucou,dthnn,dsigelm    &
     &                ,wf,wfmtp,ffrg,nrmax,i1st,fj,ls,lmax,icpot    &
     &                ,ngth,ngph,ielm,kinelm,sig                    &
     &                ,ipot,facv,facw,facvs,facws,eelm,iex,io       &
     &                ,izsp,ap,asp,aa,ab,th,costh                   &
     &                ,fnl,fnltbl,ims,iedg,ielmedg                  &
     &                ,iepfix,epfix,nepfmax,epfmin                  &
     &                ,detbl,uopttbl,ulstbl,uopt                    &
     &                ,potv,potw,potd,potrs,potis                   &
     &                ,fm0,fm1,fm2,fmsp                             &
     &                ,nthnnmax,thnnmin,nennmax,sigelm              &
     &                ,fk0,th0,ph0,fk1,th1,ph1,fk2,th2,ph2          &
     &                ,thcalg,thweig,nthgmax                        &
     &                ,phcalg,phweig,nphgmax                        &
     &                ,rcalg,rweig,nrgmax                           &
     &                ,pltbl,dpltbl,dzpl,nzplmax                    &
     &                ,hv0,c0r,fmuelm,const,tdxr                    &
     &                ,itmdcal,ilscal,wfls,uls,ay)

            sumph1a=sumph1a+tdxr*ph1weig(iph1g)*fmol

           end do

           sumk1a=sumk1a+sumph1a*fk1a*e2a*fk1weig(ik1ag)

          end do

          ddx(ikbaz,ikbab)=sumk1a/qa/hc**2*facid/fkfac**3
          if(ixunt==1) ddx(ikbaz,ikbab)=ddx(ikbaz,ikbab)*1.0d3

          sumkbb=sumkbb+ddx(ikbaz,ikbab)*fkbab*fkfac

         end do

         write(kibtbl,611) fkbazo*fkfac,(ddx(ikbaz,i),i=1,nkbabmax)
  611    format(f12.5,401e13.5)

         sumsig=sumsig+sumkbb
         pmd=sumkbb*dkbab*fkfac*2.d0*pi

         if(kiblg>0) write(kiblg,612) fkbazo*fkfac,pmd
  612    format(f12.5,2x,e13.5)

        end do

        if(kibpx+kibtr+kibtl>0) &
     &    call mdcal(nkbabmax,dkbab,nkbazmax,fkbazmin,dkbaz,kunt,ddx)

        sumsig=sumsig*dkbaz*dkbab*fkfac**2*2.d0*pi
        write(kibout,650) sumsig,csunit
  650   format(/3x,'-- integrated cross section --' &
     &         /8x,f12.2,1x,a2)
        if(nint(ap)==1 .and. nint(asp)==1) then
          call signn(elab,sigpp,sigpn)
          signntot=sigpp
          if(abs(zp-zsp)>epslim) signntot=sigpn
          write(kibout,651) signntot
  651     format(6x,'(cf. NN total cross section:',f7.2,' mb)')
        end if

       else

        etoti=eal+e0l
        fktoti=fkal+fk0l

        read(ivar,*)

        t1lprev=-999.0d0
        th1lprev=-999.0d0
        ph1lprev=-999.0d0
        t2lprev=-999.0d0
        th2lprev=-999.0d0
        ph2lprev=-999.0d0
        fkblprev=-999.0d0
        thblprev=-999.0d0
        phblprev=-999.0d0

        do
         tdx=0.0d0
         ay=0.0d0
         read(ivar,501,iostat=ios) t1l,th1l(1),ph1l(1), t2l,th2l(1),ph2l &
     &                            ,fkbl(1),thbl(1),phbl(1)
  501    format(9f11.0)
         if(kunt==1) then
           fkbl(1)=fkbl(1)/hc
          else if(kunt==2) then
           fkbl(1)=fkbl(1)/hc*1.0d3
         end if
         if(ios<0) exit
         if(ios>0) then
           write(*,*) 'ERROR: invalid format for the kinematics profile'
           stop 1
         end if

         if(ivar>40) then
           if(abs(t1l-t1lprev)<epslim .and. abs(th1l(1)-th1lprev)<epslim .and.      &
     &        abs(ph1l(1)-ph1lprev)<epslim .and. abs(t2l-t2lprev)<epslim .and.      &
     &        abs(th2l(1)-th2lprev)<epslim .and. abs(ph2l-ph2lprev)<epslim .and.    &
     &        abs(fkbl(1)-fkblprev)<epslim .and. abs(thbl(1)-thblprev)<epslim .and. &
     &        abs(phbl(1)-phblprev)<epslim) cycle
          else if(ivar>30) then
           if(abs(t1l-t1lprev)<epslim .and. abs(th1l(1)-th1lprev)<epslim .and. &
     &        abs(ph1l(1)-ph1lprev)<epslim .and. abs(t2l-t2lprev)<epslim .and. &
     &        abs(ph2l-ph2lprev)<epslim) cycle
          else if(ivar>20) then
           if(abs(th2l(1)-th2lprev)<epslim .and. abs(ph2l-ph2lprev)<epslim .and.    &
     &        abs(fkbl(1)-fkblprev)<epslim .and. abs(thbl(1)-thblprev)<epslim .and. &
     &        abs(phbl(1)-phblprev)<epslim) cycle
          else if(ivar>10) then
           if(abs(t1l-t1lprev)<epslim .and. abs(th1l(1)-th1lprev)<epslim .and. &
     &        abs(ph1l(1)-ph1lprev)<epslim .and.                               &
     &        abs(th2l(1)-th2lprev)<epslim .and. abs(ph2l-ph2lprev)<epslim) cycle
         end if

         t1lprev=t1l
         th1lprev=th1l(1)
         ph1lprev=ph1l(1)
         t2lprev=t2l
         th2lprev=th2l(1)
         ph2lprev=ph2l
         fkblprev=fkbl(1)
         thblprev=thbl(1)
         phblprev=phbl(1)

         if(ivar>=40) then
           e1l(1)=t1l+fm1
           e2l(1)=t2l+fm2

!-----     e-p conservation check

           ebl=sqrt((hc*fkbl(1))**2+fmb**2)
           fk1l=sqrt((e1l(1)+fm1)*(e1l(1)-fm1))/hc
           fk2l=sqrt((e2l(1)+fm2)*(e2l(1)-fm2))/hc
           call vecxyz(1,fk1l,   th1l(1),ph1l(1), fk1lx,fk1ly,fk1lz)
           call vecxyz(1,fk2l,   th2l(1),ph2l,    fk2lx,fk2ly,fk2lz)
           call vecxyz(1,fkbl(1),thbl(1),phbl(1), fkblx,fkbly,fkblz)

           if(abs(etoti-e1l(1)-e2l(1)-ebl)>etoti*chkcons .or. &
     &        abs(fk1lx+fk2lx+fkblx)>chkcons .or.             &
     &        abs(fk1ly+fk2ly+fkbly)>chkcons .or.             &
     &        abs(fk1lz+fk2lz+fkblz-fktoti)>fktoti*chkcons) then
             nsol=0
             if(ions==0) cycle
             call tblout(e1l(1),th1l(1),ph1l(1),e2l(1),th2l(1),ph2l    &
     &                  ,fkbl(1),thbl(1),phbl(1)                       &
     &                  ,ikin,ifrm,imir,ivar,0,0,fm1,fm2,fmb,kunt,ions &
     &                  ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
             io=io+1
             if(io*1.22d-4>fslim) then
               write(kibout,690)
               go to 1000
             end if
             cycle
           end if
           nsol=1

           iapp(1)=1
           iapp(2)=0

          else if(ivar>=30) then

           varl=t1l
           thxl=th1l(1)
           phxl=ph1l(1)
           et2l=t2l

           call kinema3bl(e0l,fk0l,eal,fkal,fmb,fm1,fm2 &
     &                   ,ivar,varl,thxl,phxl,et2l,ph2l &
     &                   ,e1l,th1l,ph1l,e2l,th2l,fkbl,thbl,phbl,nsol,iapp)

           if(nsol==0) then
             if(ions==0) cycle
             call tblout(e1l(1),th1l(1),ph1l(1),e2l(1),th2l(1),ph2l       &
     &                  ,fkbl(1),thbl(1),phbl(1)                          &
     &                  ,ikin,ifrm,imir,ivar,nsol,0,fm1,fm2,fmb,kunt,ions &
     &                  ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
             io=io+1
             if(io*1.22d-4>fslim) then
               write(kibout,690)
               go to 1000
             end if
             cycle
           end if

          else if(ivar>=20) then

           varl=fkbl(1)
           thxl=thbl(1)
           phxl=phbl(1)
           et2l=th2l(1)

           call kinema3bl(e0l,fk0l,eal,fkal,fmb,fm1,fm2 &
     &                   ,ivar,varl,thxl,phxl,et2l,ph2l &
     &                   ,e1l,th1l,ph1l,e2l,th2l,fkbl,thbl,phbl,nsol,iapp)

           if(nsol==0) then
             if(ions==0) cycle
             call tblout(e1l(1),th1l(1),ph1l(1),e2l(1),th2l(1),ph2l       &
     &                  ,fkbl(1),thbl(1),phbl(1)                          &
     &                  ,ikin,ifrm,imir,ivar,nsol,0,fm1,fm2,fmb,kunt,ions &
     &                  ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
             io=io+1
             if(io*1.22d-4>fslim) then
               write(kibout,690)
               go to 1000
             end if
             cycle
           end if

          else if(ivar>=10) then

           varl=t1l
           thxl=th1l(1)
           phxl=ph1l(1)
           et2l=th2l(1)

           call kinema3bl(e0l,fk0l,eal,fkal,fmb,fm1,fm2 &
     &                   ,ivar,varl,thxl,phxl,et2l,ph2l &
     &                   ,e1l,th1l,ph1l,e2l,th2l,fkbl,thbl,phbl,nsol,iapp)

           if(nsol==0) then
             if(ions==0) cycle
             call tblout(e1l(1),th1l(1),ph1l(1),e2l(1),th2l(1),ph2l       &
     &                  ,fkbl(1),thbl(1),phbl(1)                          &
     &                  ,ikin,ifrm,imir,ivar,nsol,0,fm1,fm2,fmb,kunt,ions &
     &                  ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)
             io=io+1
             if(io*1.22d-4>fslim) then
               write(kibout,690)
               go to 1000
             end if
             cycle
           end if

         end if

         do isol=1,2
          if(iapp(isol)==0) cycle

!-----    approved

          io=io+1

          call ltrz(1,e1l(isol), 1,th1l(isol),ph1l(isol),betgl,gamgl &
     &             ,fm1,e1,fk1,th1,ph1)
          call ltrz(1,e2l(isol), 1,th2l(isol),ph2l,betgl,gamgl       &
     &             ,fm2,e2,fk2,th2,ph2)
          call ltrz(2,fkbl(isol),1,thbl(isol),phbl(isol),betgl,gamgl &
     &             ,fmb,eb,fkb,thb,phb)

          fk1l=sqrt((e1l(isol)+fm1)*(e1l(isol)-fm1))/hc
          fk2l=sqrt((e2l(isol)+fm2)*(e2l(isol)-fm2))/hc

          if(ikin==0) then
            t1a=e1l(isol)-fm1
            t2a=e2l(isol)-fm2
           else
            call ltrz(1,e1, 1,th1,ph1,betag,gamag,fm1,e1a,fk1a,th1a,ph1a)
            call ltrz(1,e2, 1,th2,ph2,betag,gamag,fm2,e2a,fk2a,th2a,ph2a)
            t1a=e1a-fm1
            t2a=e2a-fm2
          end if

          call vecxyz(1,fkb,thb,phb, fkbx,fkby,fkbz)
          fknz=-fkbz-(aa-asp)/aa*fk0z
          fknb=sqrt(fkbx**2+fkby**2)
          fkn=sqrt(fknz**2+fknb**2)
          if(fkncut>epslim .and. fkn>fkncut) then
            if(ions==0) then
              io=io-1
              cycle
            end if
            call tblout(e1l(isol),th1l(isol),ph1l(isol)                     &
     &                 ,e2l(isol),th2l(isol),ph2l                           &
     &                 ,fkbl(isol),thbl(isol),phbl(isol)                    &
     &                 ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb,kunt,ions &
     &                 ,betgl,gamgl,betvl,gamvl,1,ical,tdx,ay)
            if(io*1.22d-4>fslim) then
              write(kibout,690)
              go to 1000
            end if
            cycle
          end if

          if(ical==0) then
            call tblout(e1l(isol),th1l(isol),ph1l(isol)                     &
     &                 ,e2l(isol),th2l(isol),ph2l                           &
     &                 ,fkbl(isol),thbl(isol),phbl(isol)                    &
     &                 ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb,kunt,ions &
     &                 ,betgl,gamgl,betvl,gamvl,0,0,tdx,ay)
            if(io*1.22d-4>fslim) then
              write(kibout,690)
              go to 1000
            end if
            cycle
          end if

          call tdxqm(rmax,dr,fkayi,thtkg,phikg                    &
     &              ,elab,amasspi,e0,e1,e2,ea,eb,t0a,t1a,t2a,ylm0 &
     &              ,zp,z1,z2,za,zb,pwsum3                        &
     &              ,rcn,ircn,pmi,beta,sfac,ucou,dthnn,dsigelm    &
     &              ,wf,wfmtp,ffrg,nrmax,i1st,fj,ls,lmax,icpot    &
     &              ,ngth,ngph,ielm,kinelm,sig                    &
     &              ,ipot,facv,facw,facvs,facws,eelm,iex,io       &
     &              ,izsp,ap,asp,aa,ab,th,costh                   &
     &              ,fnl,fnltbl,ims,iedg,ielmedg                  &
     &              ,iepfix,epfix,nepfmax,epfmin                  &
     &              ,detbl,uopttbl,ulstbl,uopt                    &
     &              ,potv,potw,potd,potrs,potis                   &
     &              ,fm0,fm1,fm2,fmsp                             &
     &              ,nthnnmax,thnnmin,nennmax,sigelm              &
     &              ,fk0,th0,ph0,fk1,th1,ph1,fk2,th2,ph2          &
     &              ,thcalg,thweig,nthgmax                        &
     &              ,phcalg,phweig,nphgmax                        &
     &              ,rcalg,rweig,nrgmax                           &
     &              ,pltbl,dpltbl,dzpl,nzplmax                    &
     &              ,hv0,c0r,fmuelm,const,tdxr                    &
     &              ,itmdcal,ilscal,wfls,uls,ay)

          phvol=0.0d0
          fjacob=1.0d0
          if(ifrm==0) then
            fktotiz=fk0l+fkal
            ebl=sqrt((hc*fkbl(isol))**2+fmb**2)
            phvol=fkin(fk1l,e1l(isol),th1l(isol),ph1l(isol),fk2l,e2l(isol)  &
     &                ,th2l(isol),ph2l,fkbl(isol),ebl,thbl(isol),phbl(isol) &
     &                ,fktotiz,ivar)
            fjacob=fjac(e1,e2,eb,e1l(isol),e2l(isol),ebl)
           else if(ifrm==1) then
            phvol=fkin(fk1,e1,th1,ph1,fk2,e2,th2,ph2,fkb,eb,thb,phb,0.0d0,ivar)
           else if(ifrm==2) then
            call ltrz(1,e1l(isol), 1,th1l(isol),ph1l(isol),betvl,gamvl &
     &               ,fm1,  e1v,fk1v,th1v,ph1v)
            call ltrz(1,e2l(isol), 1,th2l(isol),ph2l,betvl,gamvl       &
     &               ,fm2,  e2v,fk2v,th2v,ph2v)
            call ltrz(2,fkbl(isol),1,thbl(isol),phbl(isol),betvl,gamvl &
     &               ,fmb,  ebv,fkbv,thbv,phbv)
            fktotiz=-fk0v-fkav
            phvol=fkin(fk1v,e1v,th1v,ph1v,fk2v,e2v,th2v,ph2v &
     &                ,fkbv,ebv,thbv,phbv,fktotiz,ivar)
            fjacob=fjac(e1,e2,eb,e1v,e2v,ebv)
          end if

          tdx=phvol*fjacob*tdxr

          if(ixunt==1) tdx=tdx*1.0d3

          call tblout(e1l(isol),th1l(isol),ph1l(isol)                     &
     &               ,e2l(isol),th2l(isol),ph2l                           &
     &               ,fkbl(isol),thbl(isol),phbl(isol)                    &
     &               ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb,kunt,ions &
     &               ,betgl,gamgl,betvl,gamvl,0,ical,tdx,ay)

          if(io*1.22d-4>fslim) then
            write(kibout,690)
            go to 1000
          end if

         end do

        end do

      end if

 1000 continue

      return

end subroutine center
!=======================================================================



!=======================================================================
subroutine potread(ipotw,iptcl,nrmax,dr                    &
     &            ,vpotw1d,wpotw1d,vlsw1d,wlsw1d,epw,rpw   &
     &            ,vpotw,wpotw,vcouw,vlsw,wlsw,ir0cl,r0clw &
     &            ,ictrclw                                 &
     &            ,iepfixw,epfixw,epfminw,nepfmaxw,detbl   &
     &            ,uopt,uls,ucou,uopttbl,ulstbl            &
     &            ,fnlwr,fnlwi,fnlwr1d,fnlwi1d,fnl,fnltbl)
      use consts
      use dims
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 epw(nepdim),rpw(nrdim)
      integer ir0cl(3)
      real*8 vpotw(nepdim,nrdim),wpotw(nepdim,nrdim),vcouw(nrdim)
      real*8 vpotw1d(nrdim),wpotw1d(nrdim)
      real*8 vlsw(nepdim,nrdim),wlsw(nepdim,nrdim)
      real*8 vlsw1d(nrdim),wlsw1d(nrdim)
      complex*16 uopt(3,nrdim),uopttbl(3,nepfdim,nrdim)
      complex*16 uls(3,nrdim),ulstbl(3,nepfdim,nrdim)
      real*8 fnlwr(nepdim,nrdim),fnlwi(nepdim,nrdim)
      real*8 fnlwr1d(nrdim),fnlwi1d(nrdim)
      complex*16 fnl(3,nrdim),fnltbl(3,nepfdim,nrdim)
      real*8 ucou(3,nrdim)

!-----------------------------------------------------------------------

      rmax=(nrmax-1)*dr
      read(ipotw,501) nepmaxw,rpmaxw,drpw,r0clw,ictrclw
  501 format(i5,3f10.0,i5)
      if(ictrclw==0) ictrclw=1
      if(r0clw>epslim) ir0cl(iptcl)=1
      nrpmaxw=nint(rpmaxw/drpw)+1
      if(nepmaxw/=1 .and. nepmaxw<4) then
        write(*,*) 'ERROR: nepmaxw must be 1 or larger than 3'
        write(*,'(a,i0)') ' nepmaxw=',nepmaxw
        stop 1
      end if
      if(rpmaxw<rmax) then
        write(*,*) 'ERROR: rpmaxw < rmax'
        write(*,'(a,f0.4,a,f0.4)') ' rpmaxw=',rpmaxw,', rmax=',rmax
        stop 1
      end if
      if(nrpmaxw>nrdim) then
        write(*,*) 'ERROR: nrpmaxw > nrdim'
        write(*,'(a,i0,a,i0)') ' nrpmaxw=',nrpmaxw,', nrdim=',nrdim
        stop 1
      end if
      if(drpw<0.0d0) then
        write(*,*) 'ERROR: drpw must be positive'
        write(*,'(a,f0.4)') ' drpw=',drpw
        stop 1
      end if
      if(nepmaxw>nepdim) then
        write(*,*) 'ERROR: nepmaxw > nepdim'
        write(*,'(a,i0,a,i0)') ' nepmaxw=',nepmaxw,', nepdim=',nepdim
        stop 1
      end if
      epoldw=-999.0d0
      epfixw=epoldw
      do irp=1,nrpmaxw
       rpw(irp)=(irp-1)*drpw
      end do

      read(ipotw,*) epw(1)
      if(epw(1)<epslim) then
        write(*,*) 'ERROR: epw(1) must be positive'
        write(*,'(a,f0.4)') ' epw(1)=',epw(1)
        stop 1
      end if
      iflgnl=0
      do irp=1,nrpmaxw
       read(ipotw,502) vpotw(1,irp),wpotw(1,irp),vlsw(1,irp),wlsw(1,irp),vcouw(irp) &
     &                ,fnlwr(1,irp),fnlwi(1,irp)
  502  format(7e20.12)
       if(abs(fnlwr(1,irp))>epslim .or. abs(fnlwi(1,irp))>epslim) iflgnl=1
      end do
      if(iflgnl==0) then
        fnlwr(1,:)=1.0d0
        fnlwi(1,:)=0.0d0
      end if
      if(vcouw(1)>epslim) ir0cl(iptcl)=2

      do iep=2,nepmaxw
       read(ipotw,*) epw(iep)
       if(epw(iep)<epoldw) then
         write(*,*) 'ERROR: energies not in ascending order'
         stop 1
       end if
       iflgnl=0
       do irp=1,nrpmaxw
        read(ipotw,502) vpotw(iep,irp),wpotw(iep,irp),vlsw(iep,irp),wlsw(iep,irp) &
     &                 ,fnlwr(iep,irp),fnlwi(iep,irp)
        if(abs(fnlwr(iep,irp))>epslim .or. abs(fnlwi(iep,irp))>epslim) iflgnl=1
       end do
       if(iflgnl==0) then
         fnlwr(iep,:)=1.0d0
         fnlwi(iep,:)=0.0d0
       end if
      end do

      if(nepmaxw==1) then
        iepfixw=1
        epfixw=epw(1)
        do ir=1,nrmax
         rcal=(ir-1)*dr
         ircal=int(rcal/drpw+1.01)
         rbs=(ircal-1)*drpw

         if(ircal==nrpmaxw) then
           vpotfw=vpotw(1,ircal)
           wpotfw=wpotw(1,ircal)
           fnlfwr=fnlwr(1,ircal)
           fnlfwi=fnlwi(1,ircal)
          else
           vpotfw=vpotw(1,ircal)+(vpotw(1,ircal+1)-vpotw(1,ircal))*(rcal-rbs)/drpw
           wpotfw=wpotw(1,ircal)+(wpotw(1,ircal+1)-wpotw(1,ircal))*(rcal-rbs)/drpw
           fnlfwr=fnlwr(1,ircal)+(fnlwr(1,ircal+1)-fnlwr(1,ircal))*(rcal-rbs)/drpw
           fnlfwi=fnlwi(1,ircal)+(fnlwi(1,ircal+1)-fnlwi(1,ircal))*(rcal-rbs)/drpw
         end if

         uopt(iptcl,ir)=dcmplx(vpotfw,wpotfw)
         fnl(iptcl,ir)=dcmplx(fnlfwr,fnlfwi)

         if(ircal==nrpmaxw) then
           vlsfw=vlsw(1,ircal)
           wlsfw=wlsw(1,ircal)
          else
           vlsfw=vlsw(1,ircal)+(vlsw(1,ircal+1)-vlsw(1,ircal))*(rcal-rbs)/drpw
           wlsfw=wlsw(1,ircal)+(wlsw(1,ircal+1)-wlsw(1,ircal))*(rcal-rbs)/drpw
         end if
         uls(iptcl,ir)=dcmplx(vlsfw,wlsfw)

         if(ir0cl(iptcl)==2) then
           if(ircal==nrpmaxw) then
             ucou(iptcl,ir)=vcouw(ircal)
            else
             ucou(iptcl,ir)=vcouw(ircal)+(vcouw(ircal+1)-vcouw(ircal))*(rcal-rbs)/drpw
           end if
         end if

        end do
        return
      end if

      iepfixw=0
      nepfmaxw=int((epw(nepmaxw)-epw(1))/detbl+1.01)

      do iepf=1,nepfmaxw
       if(iepf==1 .and. ir0cl(iptcl)==2) then
         do ir=1,nrmax
          rcal=(ir-1)*dr
          ircal=int(rcal/drpw+1.01)
          rbs=(ircal-1)*drpw
          if(ircal==nrpmaxw) then
            ucou(iptcl,ir)=vcouw(ircal)
           else
            ucou(iptcl,ir)=vcouw(ircal)+(vcouw(ircal+1)-vcouw(ircal))*(rcal-rbs)/drpw
          end if
         end do
       end if

       epfcal=(iepf-1)*detbl+epw(1)
       epfminw=epw(1)
       do irp=1,nrpmaxw
        vpotw1d(irp)=suphodx2(epfcal,epw,vpotw,nepmaxw,irp,nepdim,nrdim)
        wpotw1d(irp)=suphodx2(epfcal,epw,wpotw,nepmaxw,irp,nepdim,nrdim)
        vlsw1d(irp) =suphodx2(epfcal,epw,vlsw, nepmaxw,irp,nepdim,nrdim)
        wlsw1d(irp) =suphodx2(epfcal,epw,wlsw, nepmaxw,irp,nepdim,nrdim)
        fnlwr1d(irp)=suphodx2(epfcal,epw,fnlwr,nepmaxw,irp,nepdim,nrdim)
        fnlwi1d(irp)=suphodx2(epfcal,epw,fnlwi,nepmaxw,irp,nepdim,nrdim)
       end do
       do ir=1,nrmax
        rcal=(ir-1)*dr
        ircal=int(rcal/drpw+1.01)
        rbs=(ircal-1)*drpw

        if(ircal==nrpmaxw) then
          vpotfw=vpotw1d(ircal)
          wpotfw=wpotw1d(ircal)
          vlsfw=vlsw1d(ircal)
          wlsfw=wlsw1d(ircal)
          fnlfwr=fnlwr1d(ircal)
          fnlfwi=fnlwi1d(ircal)
         else
          vpotfw=vpotw1d(ircal)+(vpotw1d(ircal+1)-vpotw1d(ircal))*(rcal-rbs)/drpw
          wpotfw=wpotw1d(ircal)+(wpotw1d(ircal+1)-wpotw1d(ircal))*(rcal-rbs)/drpw
          vlsfw=vlsw1d(ircal)+(vlsw1d(ircal+1)-vlsw1d(ircal))*(rcal-rbs)/drpw
          wlsfw=wlsw1d(ircal)+(wlsw1d(ircal+1)-wlsw1d(ircal))*(rcal-rbs)/drpw
          fnlfwr=fnlwr1d(ircal)+(fnlwr1d(ircal+1)-fnlwr1d(ircal))*(rcal-rbs)/drpw
          fnlfwi=fnlwi1d(ircal)+(fnlwi1d(ircal+1)-fnlwi1d(ircal))*(rcal-rbs)/drpw
        end if

        uopttbl(iptcl,iepf,ir)=dcmplx(vpotfw,wpotfw)
        ulstbl(iptcl,iepf,ir)=dcmplx(vlsfw,wlsfw)
        fnltbl(iptcl,iepf,ir)=dcmplx(fnlfwr,fnlfwi)

        if(ircal==nrpmaxw) then
          vlsfw=vlsw1d(ircal)
          wlsfw=wlsw1d(ircal)
         else
          vlsfw=vlsw1d(ircal)+(vlsw1d(ircal+1)-vlsw1d(ircal))*(rcal-rbs)/drpw
          wlsfw=wlsw1d(ircal)+(wlsw1d(ircal+1)-wlsw1d(ircal))*(rcal-rbs)/drpw
        end if
        uls(iptcl,ir)=dcmplx(vlsfw,wlsfw)

       end do
      end do

      return

end subroutine potread
!=======================================================================



!=======================================================================
subroutine tblout(e1i,th1i,ph1i,e2i,th2i,ph2i,fkbi,thbi,phbi          &
     &           ,ikin,ifrm,imir,ivar,nsol,isol,fm1,fm2,fmb,kunt,ions &
     &           ,betgl,gamgl,betvl,gamvl,iskp,ical,tdx,ay)
      use consts
      use dims
      use kibmod,only:kibtbl
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      if(kunt==0) then
        fkfac=1.0d0
       else if(kunt==1) then
        fkfac=hc
       else if(kunt==2) then
        fkfac=hc*1.0d-3
      end if
      csfac=1.0d0
      if(ivar==2) csfac=fkfac**3

      if(nsol==0) then

        if(ivar==1 .or. (ivar>=10 .and. ivar<20)) then

          if(ifrm==0) then

            if(imir==1) then
              th1i=180.0d0-th1i
              th2i=180.0d0-th2i
            end if
            write(kibtbl,601) e1i-fm1,th1i,ph1i,  th2i,ph2i,  0
  601       format(3(f10.4,1x),11x,2(f10.4,1x),44x,i3,5x,'not allowed')

           else if(ifrm==1) then

            call ltrz(1,e1i,1,th1i,ph1i,betgl,gamgl,fm1,e1o,fk1o,th1o,ph1o)
            if(imir==1) th1o=180.0d0-th1o
            write(kibtbl,602) e1o-fm1,th1o,ph1o,  0
  602       format(3(f10.4,1x),77x,i3,5x,'not allowed')

           else if(ifrm==2) then

            call ltrz(1,e1i,1,th1i,ph1i,betvl,gamvl,fm1,e1o,fk1o,th1o,ph1o)
            if(imir==1) th1o=180.0d0-th1o
            write(kibtbl,602) e1o-fm1,th1o,ph1o,  0

          end if

         else if(ivar==2 .or. (ivar>=20 .and. ivar<30)) then

          if(ikin==0) then
            fkbo=fkbi
            thbo=thbi
           else
            call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
          end if
          if(imir==1) thbo=180.0d0-thbo
          pr=fkbo*hc
          if(cos(thbo*piad)<0.0d0) pr=-pr

          if(ifrm==0) then

            if(imir==1) then
              th2i=180.0d0-th2i
              thbi=180.0d0-thbi
            end if

            write(kibtbl,603) th2i,ph2i,  fkbi*fkfac,thbi,phbi,pr,  0
  603       format(44x,6(f10.4,1x),i3,5x,'not allowed')

           else if(ifrm==1) then

            call ltrz(2,fkbi,1,thbi,phbi,betgl,gamgl,fmb,ebo,fkbo,thbo,phbo)
            if(imir==1) thbo=180.0d0-thbo
            write(kibtbl,604) fkbo*fkfac,thbo,phbo,pr,  0
  604       format(66x,4(f10.4,1x),i3,5x,'not allowed')

           else if(ifrm==2) then

            call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
            if(imir==1) thbo=180.0d0-thbo
            write(kibtbl,604) fkbo*fkfac,thbo,phbo,pr,  0

          end if

         else if(ivar==3 .or. (ivar>=30 .and. ivar<40)) then

          if(ifrm==0) then

            if(imir==1) th1i=180.0d0-th1i
            write(kibtbl,610) e1i-fm1,th1i,ph1i,e2i-fm2,  ph2i,  0
  610       format(4(f10.4,1x),11x,f10.4,1x,44x,i3,5x,'not allowed')

           else if(ifrm==1) then

            call ltrz(1,e1i,1,th1i,ph1i,betgl,gamgl,fm1,e1o,fk1o,th1o,ph1o)
            if(imir==1) th1o=180.0d0-th1o
            write(kibtbl,611) e1o-fm1,th1o,ph1o,  0
  611       format(3(f10.4,1x),77x,i3,5x,'not allowed')

           else if(ifrm==2) then

            call ltrz(1,e1i,1,th1i,ph1i,betvl,gamvl,fm1,e1o,fk1o,th1o,ph1o)
            if(imir==1) th1o=180.0d0-th1o
            write(kibtbl,602) e1o-fm1,th1o,ph1o,  0

          end if

         else if(ivar>=40) then

          if(ikin==0) then
            fkbo=fkbi
            thbo=thbi
           else
            call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
          end if
          if(imir==1) thbo=180.0d0-thbo
          pr=fkbo*hc
          if(cos(thbo*piad)<0.0d0) pr=-pr

          if(ifrm==0) then
            e1o=e1i
            th1o=th1i
            ph1o=ph1i
            e2o=e2i
            th2o=th2i
            ph2o=ph2i
            fkbo=fkbi
            thbo=thbi
            phbo=phbi
           else if(ifrm==1) then
            call ltrz(1,e1i, 1,th1i,ph1i,betgl,gamgl,fm1,e1o,fk1o,th1o,ph1o)
            call ltrz(1,e2i, 1,th2i,ph2i,betgl,gamgl,fm2,e2o,fk2o,th2o,ph2o)
            call ltrz(2,fkbi,1,thbi,phbi,betgl,gamgl,fmb,ebo,fkbo,thbo,phbo)
           else if(ifrm==2) then
            call ltrz(1,e1i, 1,th1i,ph1i,betvl,gamvl,fm1,e1o,fk1o,th1o,ph1o)
            call ltrz(1,e2i, 1,th2i,ph2i,betvl,gamvl,fm2,e2o,fk2o,th2o,ph2o)
            call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
          end if
          if(imir==1) then
            th1o=180.0d0-th1o
            th2o=180.0d0-th2o
            thbo=180.0d0-thbo
          end if
          write(kibtbl,605) e1o-fm1,th1o,ph1o,e2o-fm2,th2o,ph2o,fkbo*fkfac,thbo,phbo,pr,0
  605     format(10(f10.4,1x),i3,5x,'not allowed')

        end if

        return

      end if


      if(ikin==0) then
        fkbo=fkbi
        thbo=thbi
       else
        call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
      end if
      if(imir==1) thbo=180.0d0-thbo
      pr=fkbo*hc
      if(cos(thbo*piad)<0.0d0) pr=-pr

      if(ifrm==0) then
        e1o=e1i
        th1o=th1i
        ph1o=ph1i
        e2o=e2i
        th2o=th2i
        ph2o=ph2i
        fkbo=fkbi
        thbo=thbi
        phbo=phbi
       else if(ifrm==1) then
        call ltrz(1,e1i, 1,th1i,ph1i,betgl,gamgl,fm1,e1o,fk1o,th1o,ph1o)
        call ltrz(1,e2i, 1,th2i,ph2i,betgl,gamgl,fm2,e2o,fk2o,th2o,ph2o)
        call ltrz(2,fkbi,1,thbi,phbi,betgl,gamgl,fmb,ebo,fkbo,thbo,phbo)
       else if(ifrm==2) then
        call ltrz(1,e1i, 1,th1i,ph1i,betvl,gamvl,fm1,e1o,fk1o,th1o,ph1o)
        call ltrz(1,e2i, 1,th2i,ph2i,betvl,gamvl,fm2,e2o,fk2o,th2o,ph2o)
        call ltrz(2,fkbi,1,thbi,phbi,betvl,gamvl,fmb,ebo,fkbo,thbo,phbo)
      end if

      if(imir==1) then
        th1o=180.0d0-th1o
        th2o=180.0d0-th2o
        thbo=180.0d0-thbo
      end if

      if(iskp==1) then
        write(kibtbl,606) &
     &    e1o-fm1,th1o,ph1o,e2o-fm2,th2o,ph2o,fkbo*fkfac,thbo,phbo,pr,-isol
  606   format(10(f10.4,1x),i3,5x,'skipped')
       else if(ical==0) then
        write(kibtbl,607) &
     &    e1o-fm1,th1o,ph1o,e2o-fm2,th2o,ph2o,fkbo*fkfac,thbo,phbo,pr,isol
  607   format(10(f10.4,1x),i3)
       else if(abs(tdx)<1.0d-99) then
        if(ions==0) return
        write(kibtbl,608) &
     &    e1o-fm1,th1o,ph1o,e2o-fm2,th2o,ph2o,fkbo*fkfac,thbo,phbo,pr,-isol
  608   format(10(f10.4,1x),i3,5x,'zero')
       else
        write(kibtbl,609) &
     &    e1o-fm1,th1o,ph1o,e2o-fm2,th2o,ph2o,fkbo*fkfac,thbo,phbo,pr,isol,tdx/csfac,ay
  609   format(10(f10.4,1x),i3,3x,e13.5,4x,f10.4)
      end if

      return

end subroutine tblout
!=======================================================================



!=======================================================================
subroutine sigelm0(ielm,elab,izp,izsp,ap,asp,nthnnmax,thnnmin,dthnn,nennmax &
     &            ,sigelm,eelm,dsigelm)
      use consts
      use dims
      use kibmod
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 eelm(nenndim),dsigelm(nenndim,nthnndim)

!-----------------------------------------------------------------------

      if(ielm==0) then

!-----  isotropic free nn total x-section at elab

        call signn(elab,sigpp,sigpn)

        if(izp==izsp) then
!-----    the pp total cross section (experimental value) is defined
!-----    by the pp differential cross section integrated over the
!-----    c.m. scattering angle from 0 to pi/2 (not to pi).
          sigelm=sigpp*2.0d0/(4.0d0*pi)/10.0d0
         else
          sigelm=sigpn/(4.0d0*pi)/10.0d0
        end if
          fmuelm=reduced_energy(elab,ap,asp)
          sigelm=sigelm*(2.0d0*pi*hc**2)**2/fmuelm**2.d0

       else if(ielm==3) then

!-----  free differential x-section (making a table)

        thnnmin=0.0d0
        eold=-999.0d0
        read(kibelm,*) nennmax,thnnmax,dthnn
        if(thnnmax<0.0d0) then
          write(*,*) 'ERROR: thnnmax must be positive'
          write(*,'(a,f0.4)') ' thnnmax=',thnnmax
          stop 1
        end if
        if(dthnn<0.0d0) then
          write(*,*) 'ERROR: dthnn must be positive'
          write(*,'(a,f0.4)') ' dthnn=',dthnn
          stop 1
        end if
        nthnnmax=nint(thnnmax/dthnn)+1

        if(nthnnmax>nthnndim) then
          write(*,*) 'ERROR: nthnnmax > nthnndim'
          write(*,'(a,i0,a,i0)') ' nthnnmax=',nthnnmax,', nthnndim=',nthnndim
          stop 1
        end if

        if(nennmax>nenndim) then
          write(*,*) 'ERROR: nennmax > nenndim'
          write(*,'(a,i0,a,i0)') ' nennmax=',nennmax,', nenndim=',nenndim
          stop 1
        end if

        ilinetot=nennmax*nthnnmax+nennmax
        iline=0
        ie=0
        do
         read(kibelm,*,iostat=ios) eelmw
         if(ios<0) goto 999
         ie=ie+1
         eelm(ie)=eelmw
         iline=iline+1
         if(eelm(ie)<eold) then
           write(*,*) 'energies not in ascending order'
           stop 1
         end if
         eold=eelm(ie)
         fmuelm=reduced_energy(eelm(ie),ap,asp)
         fmassfac=(2.0d0*pi*hc**2.d0)**2.d0/fmuelm**2.d0
         do ith=1,nthnnmax
          read(kibelm,*,end=999) dsigpp,dsigpn
          iline=iline+1
          dsigelm(ie,ith)=dsigpp/10.0d0 *fmassfac
          if(nint(ap+asp)==2 .and. izp+izsp==1) then
            dsigelm(ie,ith)=dsigpn/10.0d0 *fmassfac
          end if
         end do

        end do
    999 continue
        if(ilinetot/=iline) then
          write(*,'(a,i0,a)') 'total line numbers of input # ',kibelm, &
     &                        ' is inconsistent with header information'
          write(*,'(a,i0)') 'lines expected = ',ilinetot+1
          write(*,'(a,i0)') 'actual lines read-in = ',iline
          stop 1
        end if

      end if

      return

end subroutine sigelm0
!=======================================================================



!=======================================================================
real*8 function reduced_energy(ekinp,ap,asp)
      use consts,only:hc,ac
      implicit real*8 (a-h,o-z)
      real*8, intent(in) :: ekinp,ap,asp
!=======================================================================

!-----------------------------------------------------------------------

      eptotlab=ekinp+ap*ac
      esptotlab=asp*ac
      fkplab=sqrt(eptotlab**2.d0-(ap*ac)**2.d0)
      fksplab=0.d0
      beta=(fkplab+fksplab)/(eptotlab+esptotlab)
      gam=1.d0/sqrt(1.d0-beta**2.d0)
      call ltrz(1,eptotlab,1,0.d0,0.d0,beta,gam,ap*ac,eptotcm,fkpcm,thpcm,phpcm)
      call ltrz(1,esptotlab,1,0.d0,0.d0,beta,gam,asp*ac,esptotcm,fkspcm,thspcm,phspcm)

      reduced_energy=(eptotcm*esptotcm)/(eptotcm+esptotcm)

end function reduced_energy
!=======================================================================



!=======================================================================
!****
!****                      ***** ylmph0 *****
!****
!****    spherical harmonics y_l^m(z) for 0 =< l =< lmax, 0 =< m =< mmax
!****             with z=cos(theta) and phi=0.0 are calculated
!****
!****                                                  2007/11/12 K.Ogata
!****
!****      ldim   : dimension of l => lmax+1
!****      mdim   : dimension of m => mmax+1
!****      nthdim : dimension of theta => nthmax+1
!****      lmax   : maximum value of l
!****      mmax   : maximum value of m
!****      nthmin : minimum value of nth
!****      nthmax : maximum value of nth
!****      costh  : array of cos(theta(nth))
!****               ylm's for costh(nthmin), costh(nthmin+1), ...,
!****               costh(nthmax) are calculated
!****      ylm0   : spherical harmonics y_l^m(z) at phi=0.0 (m >= 0)
!****
subroutine ylmph0(ldim,mdim,nthdim,lmax,mmax,nthmin,nthmax,costh,ylm0)
      use consts
      use kibmod
      use angmom,only:faclog,dfaclog,memoam
      implicit real*8(a-h,o-z)
!=======================================================================
      real*8 ylm0(nthdim,ldim,mdim),costh(nthdim)

!-----------------------------------------------------------------------

      lmax1=lmax+1
      mmax1=mmax+1

      if(lmax1<=ldim.and.mmax1<=mdim) go to 100
      write(kibout,101) lmax,ldim,mmax,mdim
  101 format(/1x,'******* dimension over in "ylmph0" --- ' &
     &          ,'lmax,ldim=',2i5,'   mmax,mdim=',2i5)
      stop 1
  100 continue


      if(memoam >= lmax+mmax+1) go to 1

!--- faclog(n) means log [(n-1)!]
      faclog(1)=0.0d0
      faclog(2)=0.0d0
      do i=3,lmax+mmax+1
       fi=i-1
       faclog(i)=faclog(i-1)+dlog(fi)
      end do

!--- dfaclog(n) means log [(2*(n-1)-1)!!]
      dfaclog(1)=0.0d0
      fn=0.0d0
      do i=2,mmax+1
       fn=fn+1.0d0
       dfaclog(i)=dfaclog(i-1)+dlog(2.0d0*fn-1.0d0)
      end do

      memoam=lmax+mmax+1

    1 continue

! -----
      do it=nthmin,nthmax

       ylm0(it,:,:)=0.0d0

       z=costh(it)

!---- starting values p_m^m

       do ll=1,lmax1
        l=ll-1
        m=l
        fl=l*1.0d0
        fm=fl

!-----  p_lm (l=m, m+1, ... l_max) is renormalized to 1/(2m-1)!! x (1 - z^2)^(m/2)

        ylm0(it,ll,ll)=(-1.0d0)**fm
        if(ll/=lmax1) ylm0(it,ll+1,ll)=z*(2.0d0*fm+1.0d0)*ylm0(it,ll,ll)

       end do

       do mm=1,lmax1
        if(mm+2>lmax1) cycle
        do ll=mm+2,lmax1
         l=ll-1
         m=mm-1
         fl=l*1.0d0
         fm=m*1.0d0

         ylm0(it,ll,mm)=( z*(2.0d0*fl-1.0d0)*ylm0(it,ll-1,mm) &
     &                    -(fl+fm-1.0d0)*ylm0(it,ll-2,mm) )/(fl-fm)

        end do
       end do

      end do

! -----

      facpi=1.0d0/dsqrt(4.0d0*pi)
      do ll=1,lmax1
       ll1=ll-1
       facl=dsqrt(2.0d0*ll1+1.0d0)
       do mm=1,mmax1
        mm1=mm-1
        facfac=elemfclog(ll1,mm1,ldim,mdim,faclog)

        do it=nthmin,nthmax
         z=costh(it)
         if(mm1==0) then
           facylmlg=facfac+dfaclog(mm)
           facylm=ylm0(it,ll,mm)*dexp(facylmlg)
          else if(abs(z**2-1.0d0)<1.0d-15) then
           facylm=0.0d0
          else
           facylmlg=facfac+dfaclog(mm)+mm1/2.0d0*dlog(1.0d0-z**2)
           facylm=ylm0(it,ll,mm)*dexp(facylmlg)
         end if
         ylm0(it,ll,mm)=facpi*facl*facylm
        end do
       end do
      end do

      return

end subroutine ylmph0
!=======================================================================



!=======================================================================
!**
!**       ******** elemfclog ********
!**             elemfclog= dlog{dsqrt( (l-m)!/(l+m)! )}
!**
function elemfclog(l,m,ldim,mdim,faclog)
      implicit  real*8(a-h,o-z)
!=======================================================================
      real*8 faclog(ldim+mdim+1)

!-----------------------------------------------------------------------

      elemfclog=0.0d0
      k1=l-m
      k2=l+m
      if(k1<0) return
      a=(faclog(k1+1)-faclog(k2+1)) /2.0d0
      elemfclog=a
      return

end function elemfclog
!=======================================================================



!=======================================================================
!**
!**
!**        ***** bound2 *****
!**
!            original program is  H.Yosida's "bound"
!                           and "ffsub4" in igarashi's "twofnr".

!      2001.3  so-pot. corrected
!      2004.8  interface modified rather drastically

subroutine bound2(ffmm,zze2,ls,nrmax,dr,energy,ffr,vtot &
     &           ,vcen,vls,vcou,radi,radc,pmass         &
     &           ,pvc,pvs,ish,nod,vdeptho,wlso)

      use consts
      use dims
      use kibmod
      implicit real*8(a-h,o-z)
      real*8 vtot(nrdim),vcen(nrdim),vls(nrdim),vcou(nrdim)
      real*8 ffr(nrdim),ast(20)

!----- ffmm = 2*mu*ac/hc**2
!----- zze2  = zp*zt * hc/ec

!** =======================================================

      vdeptho=0.0d0
      wlso=0.0d0

      fll=ls*1.01d0

      fcls=0.d0
      if(ish==0) then
        ktlvbe=2
        vdepth=1.0d0
        wls=1.0d0
       else if(ish==1) then
        ktlvbe=1
        bengy=-pvc
        wls=1.0d0
       else if(ish==2) then
        ktlvbe=3
        bengy=-pvc
        fcls =-pvs
       else
        write(kibout,*) 'error in ish'
        stop 1
      end if

      eps7=1.0d-8
      test=1.0d16
      lmom=nint(fll)
      rmax=dr*(nrmax-1)
      lmom1=lmom+1
      lmom2=lmom1+1
      flmom=lmom
      flmom1=lmom1
      fl1=flmom*flmom1

      radz=radi
      if(radz==0.0d0) radz=1.2d0

      do 2000 ii=2,30

      drz= 0.3d0*dble(ii/2)*(-1.0)**ii

      korec=0
      niter=0
      incr=0
      if(ktlvbe-2) 10,20,10
   10  vdepth=bengy +(3.1415926d0*( nod+flmom1/2))**2 &
     &                /(0.048228d0*pmass*(radi+drz)**2)
       if(ktlvbe==3) wls=fcls*vdepth
       go to 80

   20  bengy=vdepth -(3.1415926d0*( nod+flmom1/2))**2 &
     &                 /(0.048228d0*pmass*radz*radz)
       if(bengy-eps7) 25,25,80
   25  if(dabs(drz)<0.01) drz=0.3
       radz=radz+drz
       incr=incr+1
       if(incr<=20) go to 20
       kcheck=11
       go to 2000

   80 continue

      wk=dsqrt(ffmm*bengy)

      wrhoc=wk*radc
      wrhoz=wk*radz
      wrhocs=wrhoc*wrhoc

      weta= zze2/2.0d0 *wk/bengy

      wetac=weta/wrhoc
      wdrho=dr*wk

      rhoa=rmax*wk
      drhosq=wdrho*wdrho
      dr56=5.0d0*drhosq/6.0d0
      dr12=0.10d0*dr56

!    search classical turning point
  100 continue
      wvc=vdepth/bengy
      wvs=wls/bengy
      i=nrmax
      x1=rhoa
      t1=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
      do 110 i=nrmax-1,2,-1
       x1=x1-wdrho
       t2=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
       if(t1*t2<=0.0) go to 115
       t1=t2
  110 continue
  115 match=i+1

      if(match<=3) match=4

!    inner solution
      zer=1.0d0
      ffr(1)=0.0d0
      wrho=wdrho
      do 180 j=2,lmom2+1
       a1=-wvs*vls(j)*wrho/(flmom1+flmom1)
       b1=1.0d0-wvc*vcen(j)+3.0d0*wetac
       b2=wvs*vls(j)*wrho
       a2=(b1-b2*a1)/(4.0d0*flmom1+2.0d0)
       a3=(b1*a1-b2*a2)/(6.0d0*flmom1+6.0d0)
       wrhosq=wrho*wrho
       b3=weta/(wrhoc*wrhocs)
       a4=(b1*a2-b2*a3-b3)/(8.0d0*flmom1+12.0d0)
       a5=(b1*a3-b2*a4-b3*a1)/(10.0d0*flmom1+20.0d0)
       a6=(b1*a4-b2*a5-b3*a2)/(12.0d0*flmom1+30.0d0)
       ffr(j)=(wrho**lmom1)*(1.0d0 +a1*wrho +a2*wrhosq +a3*wrho*wrhosq &
     &        +a4*wrhosq*wrhosq +a5*wrho*wrhosq*wrhosq +a6*wrhosq**3)
       wrho=wrho+wdrho
  180 continue

      i=lmom2
      x1=wdrho*(i-1)
      t1=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
      i=i+1
      x1=x1+wdrho
      t2=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
      do 200 i=lmom2+2,match+3
       x1=x1+wdrho
       t3=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
       fac1= 1.0d0 -dr12*t1
       fac2= 2.0d0 +dr56*t2
       fac3= 1.0d0 -dr12*t3
       ffr(i)=(ffr(i-1)*fac2-ffr(i-2)*fac1)/fac3
       t1=t2
       t2=t3
       if(dabs(ffr(i))>test) then
         do 190 ip=1,i
           ffr(ip)=ffr(ip)/test
  190    continue
         zer=zer/test
       end if
  200 continue

!    node check
      nodes= 0
      do 220 i=2,match
       if(ffr(i)*ffr(i+1)) 210,215,220
  210  nodes=nodes+2
       go to 220
  215  nodes=nodes+1
  220 continue
      nn=nodes/2
      if(nn==nod) go to 240

      korec=korec+1
      if(korec>150) then
        kcheck=10
        go to 2000
      end if
      vcor= (wrhoz*wrhoz+9.86959d0*(nod+0.5d0*flmom1)**2) &
     &    / (wrhoz*wrhoz+9.86959d0*(nn +0.5d0*flmom1)**2)
      if(vcor>1.2) vcor=1.20
      if(vcor<0.8) vcor=0.80

      if(ktlvbe/=2) then
        vdepth=vcor*vdepth
        if(ktlvbe==3) wls=vcor*wls
        go to 100
      else
        bengy=bengy/vcor
        go to 80
      end if

  240 ffim3=ffr(match-3)
      ffim2=ffr(match-2)
      ffim1=ffr(match-1)
      ffi0 =ffr(match)
      ffip1=ffr(match+1)
      ffip2=ffr(match+2)
      ffip3=ffr(match+3)
      dffr1=((ffip3-ffim3)/60.0d0 -3.0d0*(ffip2-ffim2)/20.0d0 &
     &          +3.0d0*(ffip1-ffim1)/4.0d0 ) / wdrho

!    outer solution
      ast(1)=1.0d0
      t1=1.0d0
      t2=2.0d0
      do 320 i=1,18
        ast(i+1)=(flmom1-t1-weta)*(flmom1+t1-1.0d0+weta)*ast(i)/t2
        t1=t1+1.0d0
        t2=t2+2.0d0
  320 continue

      wrho=rhoa-wdrho
      do 340 j=nrmax-1,nrmax
      frt=1.0d0
      trhoa=wrho
      do 330 i=1,18
        frt=frt+ast(i+1)/trhoa
        trhoa=trhoa*wrho
  330 continue

      ffr(j)=frt/dexp(wrho+weta*dlog(wrho+wrho))
      wrho=wrho+wdrho
  340 continue

      i=nrmax
      x1=rhoa
      t1=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
      i=i-1
      x1=x1-wdrho
      t2=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
      do 360 i=nrmax-2,match-3,-1
       x1=x1-wdrho
       t3=1.0d0-wvc*vcen(i)-wvs*vls(i)+vcou(i)/bengy+fl1/(x1*x1)
       fac1=1.0d0-dr12*t1
       fac2=2.0d0+dr56*t2
       fac3=1.0d0-dr12*t3
       ffr(i)=(ffr(i+1)*fac2-ffr(i+2)*fac1)/fac3
       t1=t2
       t2=t3
       if(dabs(ffr(i))>test) then
        do 356 ip=nrmax,i,-1
         ffr(ip)=ffr(ip)/test
  356   continue
       end if
  360 continue

      ffom3=ffr(match-3)
      ffom2=ffr(match-2)
      ffom1=ffr(match-1)
      ffo0 =ffr(match)
      ffop1=ffr(match+1)
      ffop2=ffr(match+2)
      ffop3=ffr(match+3)
      dffr2=((ffop3-ffom3)/60.0d0 -3.0d0*(ffop2-ffom2)/20.0d0 &
     &          +3.0d0*(ffop1-ffom1)/4.0d0 ) / wdrho

!    matching
      ffr(match-3)=ffim3
      ffr(match-2)=ffim2
      ffr(match-1)=ffim1
      ffr(match)  =ffi0

      ratio= ffi0/ffo0
      tlogd1=dffr1/ffi0
      tlogd2=dffr2/ffo0
      difnce=dabs(tlogd1-tlogd2)
      if(difnce<=eps7) go to 510
      niter=niter+1
      if(niter>60) then
         kcheck=12
         go to 2000
      end if

      ratio2=ratio*ratio
      fnum=ffo0*dffr2*ratio2 - ffi0*dffr1
      sum=0.0d0
      if(ktlvbe/=2) then
        do 420 i=1,match
          sum= sum +ffr(i)*ffr(i)*vcen(i)*wvc
  420   continue
        do 421 i=match+1,nrmax
          sum= sum +ffr(i)*ffr(i)*vcen(i)*wvc*ratio2
  421   continue
        sum=-sum
      else
        do 430 i=1,match
          sum= sum +ffr(i)*ffr(i)
  430   continue
        do 431 i=match+1,nrmax
          sum= sum +ffr(i)*ffr(i)*ratio2
  431   continue
      end if

      denom=sum*wdrho
      incr=0
      ram1=fnum/denom
  482 ramda=1.0d0+ram1
      if(ramda-eps7) 485,485,488
  485 ram1=0.5*ram1
      incr=incr+1
      if(incr<=20) go to 482
        kcheck=13
        go to 2000
  488 korec=0
      if(ramda>1.1) ramda=1.1
      if(ramda<0.9) ramda=0.9
      if(ktlvbe/=2) then
        vdepth=ramda*vdepth
        if(ktlvbe==3) wls=ramda*wls
        go to 100
       else
        bengy=bengy*ramda
        go to 80
      end if

 2000 continue
      go to 7400

!    solution is found

  510 kcheck=0
      do 520 i=match+1,nrmax
        ffr(i)=ffr(i)*ratio
  520 continue

      sum=0.0d0
      do 570 i=1,nrmax
        sum=sum+ffr(i)*ffr(i)
  570 continue
      sum=sum*dr
      znorm=1.0d0/dsqrt(sum)
      r=0.0d0
      do 600 i=2,nrmax
        r=r+dr
        ffr(i)=ffr(i)*znorm/r
  600 continue

!-----for l=0, the radial w.fn. goes to [sin(kr)]/r \sim k around r=0.
!-----the overall magnitude is determined by znorm.
      if(lmom==0) ffr(1)=znorm*wk
      if(ffr(2)<0.0d0) ffr(1)=-ffr(1)

      if(ktlvbe/=2) then
       do i=1,nrmax
        vcen(i)=vdepth*vcen(i)
        if(ktlvbe==3) vls(i)=wls*vls(i)
        vtot(i)=-vcen(i)-vls(i)+vcou(i)
       end do
      end if
      energy=-bengy
      ish=0


!*****output*****

      if(ktlvbe==1) then
        vdeptho=vcen(1)
        wlso=pvs
       else if(ktlvbe==2) then
        vdeptho=pvc
        wlso=pvs
       else
        vdeptho=vdepth
        wlso=wls
      end if

!     ratio1=wlso/vdeptho

      write(kibout,681) bengy,vdeptho
  681 format(/3x,'-- bound-state calculation outputs --' &
     &       /8x,'binding energy               :',f10.5  &
     &       /8x,'central potential depth      :',f10.5)

      if(abs(wlso)<epslim) return
      if(lmom/=0) then
        write(kibout,682) wlso
  682   format(8x,'spin-orbit potential depth   :',f10.5)
       else
        write(kibout,683) wlso
  683   format(8x,'spin-orbit potential depth   :',f10.5 &
     &           ,' (but not-effective)')
      end if
      return


 7400 write(kibout,7410) &
     &   kcheck,korec,niter,ramda,vcor,tlogd1,tlogd2,ratio
 7410 format(/1x ,'***no convergence in "bound" -- kcheck=',i3 &
     &       /1x ,10x,'korec, niter =',2i4                     &
     &       /1x ,10x,'ramda, vcor  =',2f10.5                  &
     &       /1x ,10x,'tlogd1, tlogd2, ratio =',1p3d13.5)
      stop 1

end subroutine bound2
!=======================================================================



!=======================================================================
subroutine potbs(izsp,asp,za,aa,fj,ls,ebind,zzeb,nrmax,dr           &
     &          ,ibmc,rc,ictrc,a0c,rcl,ictrcl,ibms,v0ls,rs,ictrs,as &
     &          ,radi,rads,radc,vtot,vcen,vls,vcou,vdepc,vdeps)
      use consts
      use dims
      implicit real*8 (a-h,o-z)

!-----potential parameter given by Bohr and Mottelson

!----- the single-particle potential of Bohr and Mottelson has
!----- the following form:
!-----   u = vc f(r) + vs (\ell \dot s) r_0^2/r df(r)/dr
!----- with
!-----   f(r) = (1 + exp[(r-r_0 A^(1/3))/a])^(-1),
!-----   r_0 = 1.27 fm,  a = 0.67 fm,
!-----   vc = -51.0 - 33.0 (N-Z)/A  MeV    for proton,
!-----   vc = -51.0 + 33.0 (N-Z)/A  MeV    for neutron,
!-----   vs = -0.44 vc.
!----- in bound2 the potential is defined by
!-----   u' = vc' f(r) + vs' 2 (\ell \dot s) 2/r
!-----                   x 1/a exp[(r-r_0 A^(1/3))/a]
!-----                         /(1 + exp[(r-r_0 A^(1/3))/a])^2.
!----- since
!-----   u = vc f(r) + vs (\ell \dot s) r_0^2/r df(r)/dr
!-----     = vc f(r) + (0.44 r_0^2/4 vc) 2 (\ell \dot s) 2/r
!-----                 x 1/a exp[(r-r_0 A^(1/3))/a]
!-----                       /(1 + exp[(r-r_0 A^(1/3))/a])^2,
!----- one finds
!-----   vc' = vc
!-----   vs' = 0.44 r_0^2/4 vc
!----- and the ratio of vs' to vc' is 0.44 r_0^2/4 (= 0.177419).
!=======================================================================
      real*8 vtot(nrdim),vcen(nrdim),vls(nrdim),vcou(nrdim)

!-----------------------------------------------------------------------

      iza=nint(za)
      iaa=nint(aa)
      ina=iaa-iza

      vtot(:)= 0.0d0
      vcen(:)= 0.0d0
      vls(:) = 0.0d0
      vcou(:)= 0.0d0

      if(ibmc==1 .or. ibms==1) then
        if(izsp==1) then
          vdepc=(-51.0d0-33.0d0*(ina-iza)*1.0d0/iaa)*(-1.0d0)
         else
          vdepc=(-51.0d0+33.0d0*(ina-iza)*1.0d0/iaa)*(-1.0d0)
        end if
        vdeps=vdepc*0.177419d0
      end if

      if(ibmc==1) then
        rc=1.27d0
        ictrc=0
        a0c=0.67d0
        rcl=rc
        ictrcl=0
       else
        vdepc=ebind
      end if
      if(ibms==1) then
        rs=1.27d0
        ictrs=0
        as=0.67d0
       else
        vdeps=v0ls
      end if

      if(ictrc==0) then
        facc=aa**(1.0d0/3.0d0)
       else if(ictrc==1) then
        facc=(aa-asp)**(1.0d0/3.0d0)
       else if(ictrc==2) then
        facc=(aa-asp)**(1.0d0/3.0d0)+asp**(1.0d0/3.0d0)
       else if(ictrc==3) then
        facc=1.0d0
      end if

      if(ictrcl==0) then
        faccl=aa**(1.0d0/3.0d0)
       else if(ictrcl==1) then
        faccl=(aa-asp)**(1.0d0/3.0d0)
       else if(ictrcl==2) then
        faccl=(aa-asp)**(1.0d0/3.0d0)+asp**(1.0d0/3.0d0)
       else if(ictrcl==3) then
        faccl=1.0d0
      end if

      if(ictrs==0) then
        facs=aa**(1.0d0/3.0d0)
       else if(ictrs==1) then
        facs=(aa-asp)**(1.0d0/3.0d0)
       else if(ictrs==2) then
        facs=(aa-asp)**(1.0d0/3.0d0)+asp**(1.0d0/3.0d0)
       else if(ictrs==3) then
        facs=1.0d0
      end if

      r0c =rc*facc
      r0cl=rcl*faccl
      r0s =rs*facs

      radi=r0c
      rads=r0s
      radc=r0cl

!*****central*****

      do nr=1,nrmax
       r=dr*(nr-1)
       x= (r-r0c)/a0c
       vcen(nr)=vcen(nr) + vdepc/(1.0d0+dexp(x))
       vtot(nr)=vtot(nr) - vdepc/(1.0d0+dexp(x))
      end do

!*****spin-orbit*****

      factls=(fj*(fj+1.0d0) -ls*(ls+1) -0.75d0)/2
      do nr=1,nrmax
       r=dr*(nr-1)
       rdv=r
       if(r<1.0d-10) rdv=1.0d-10
       x= (r-r0s)/as
       vls(nr)=vls(nr)+2.0d0*factls/0.5d0/rdv*vdeps/as*dexp(x)/(1.0d0+dexp(x))**2
       vtot(nr)=vtot(nr)-2.0d0*factls/0.5d0/rdv*vdeps/as*dexp(x)/(1.0d0+dexp(x))**2
      end do

!*****coulomb*****

      do nr=1,nrmax
       r=dr*(nr-1)
       if(nr==1) r=1.0d-10
       vcou(nr) =zzeb/r
       if(r<radc) vcou(nr) =zzeb/2.0*(3.0-(r/radc)**2)/radc
       vtot(nr)=vtot(nr) + vcou(nr)
      end do

      return

end subroutine potbs
!=======================================================================



!=======================================================================
subroutine signn(elab,sigpp,sigpn)
      implicit real*8 (a-h,o-z)

!-----parameterized by C.A. Bertulani and C.De Conti
!-----Phys. Rev. C 81, 064603 (2010)
!=======================================================================

!-----------------------------------------------------------------------

      if(elab<=0.0d0 .or. elab>5.0d3) then
        write(*,*) 'ERROR: elab is out of range'
        write(*,'(a,e13.5)') ' elab=',elab
        stop 1
      end if

      if(elab<280.0d0) then
        sigpp=19.6d0+4253.0d0/elab-375.0d0/sqrt(elab)+3.86d-2*elab
       else if(elab<840.0d0) then
        sigpp=32.7d0-5.52d-2*elab+3.53d-7*elab**3-2.97d-10*elab**4
       else
        sigpp=50.9d0-3.8d-3*elab+2.78d-7*elab**2+1.92d-15*elab**4
      end if

      if(elab<300.0d0) then
        sigpn=89.4d0-2025.0d0/sqrt(elab)+19108.0d0/elab-43535.0d0/elab**2
       else if(elab<700.0d0) then
        sigpn=14.2d0+5436.0d0/elab+3.72d-5*elab**2-7.55d-9*elab**3
       else
        sigpn=33.9d0+6.1d-3*elab-1.55d-6*elab**2+1.3d-10*elab**3
      end if

      return

end subroutine signn
!=======================================================================



!=======================================================================
subroutine kinema3bl(e0l,fk0l,eal,fkal,fmb,fm1,fm2 &
     &              ,ivar,varl,thxl,phxl,et2l,ph2l &
     &              ,e1l,th1l,ph1l,e2l,th2l,fkbl,thbl,phbl,nsol,iapp)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 e1l(2),th1l(2),ph1l(2),e2l(2),th2l(2),fkbl(2),thbl(2),phbl(2)
      integer iapp(2)

!-----------------------------------------------------------------------

      nsol=0
      iapp(1)=0
      iapp(2)=0

      if(ivar==1 .or. (ivar>=10 .and. ivar<20)) then
        e1l(1)=varl+fm1
        th1l(1)=thxl
        ph1l(1)=phxl
        e1l(2)=e1l(1)
        th2l(1)=et2l
        th1l(2)=th1l(1)
        ph1l(2)=ph1l(1)
        th2l(2)=th2l(1)

        fk1l=sqrt((e1l(1)+fm1)*(e1l(1)-fm1))/hc

        call vecxyz(1,fk1l,th1l(1),ph1l(1), fk1lx,fk1ly,fk1lz)
        call vecxyz(1,1.0d0,th2l(1),ph2l, uk2lx,uk2ly,uk2lz)

        qx=-fk1lx
        qy=-fk1ly
        qz=fk0l+fkal-fk1lz
        q2=qx**2+qy**2+qz**2

        ei1=e0l+eal-e1l(1)
        faca=0.5d0*(ei1**2+fm2**2-fmb**2-hc**2*q2)
        facb=hc**2*(qx*uk2lx+qy*uk2ly+qz*uk2lz)
        facc=facb**2-ei1**2*hc**2

        d2=faca**2*facb**2-(faca**2-ei1**2*fm2**2)*facc
        if(d2<0.0d0) return

        do isol=1,2

         if(isol==1) then
           fk2l=(-faca*facb-sqrt(d2))/facc
           if(fk2l<epslim) cycle
          else
           fk2l=(-faca*facb+sqrt(d2))/facc
           if(fk2l<epslim) cycle
         end if

         iapp(isol)=1
         e2l(isol)=sqrt((hc*fk2l)**2+fm2**2)

         fkblx=qx-fk2l*uk2lx
         fkbly=qy-fk2l*uk2ly
         fkblz=qz-fk2l*uk2lz

         call vecthph(1,fkblx,fkbly,fkblz,fkbl(isol),thbl(isol),phbl(isol))

        end do

       else if(ivar==2 .or. (ivar>=20 .and. ivar<30)) then

        fkbl(1)=varl
        thbl(1)=thxl
        phbl(1)=phxl
        th2l(1)=et2l
        fkbl(2)=fkbl(1)
        thbl(2)=thbl(1)
        phbl(2)=phbl(1)
        th2l(2)=th2l(1)

        ebl=sqrt((hc*fkbl(1))**2+fmb**2)

        call vecxyz(1,fkbl(1),thbl(1),phbl(1), fkblx,fkbly,fkblz)
        call vecxyz(1,1.0d0,th2l(1),ph2l, uk2lx,uk2ly,uk2lz)

        qx=-fkblx
        qy=-fkbly
        qz=fk0l+fkal-fkblz
        q2=qx**2+qy**2+qz**2

        eib=e0l+eal-ebl
        faca=0.5d0*(eib**2+fm2**2-fm1**2-hc**2*q2)
        facb=hc**2*(qx*uk2lx+qy*uk2ly+qz*uk2lz)
        facc=facb**2-eib**2*hc**2

        d2=faca**2*facb**2-(faca**2-eib**2*fm2**2)*facc
        if(d2<0.0d0) return

        do isol=1,2

         if(isol==1) then
           fk2l=(-faca*facb-sqrt(d2))/facc
           if(fk2l<epslim) cycle
          else
           fk2l=(-faca*facb+sqrt(d2))/facc
           if(fk2l<epslim) cycle
         end if

         iapp(isol)=1
         e2l(isol)=sqrt((hc*fk2l)**2+fm2**2)

         fk1lx=qx-fk2l*uk2lx
         fk1ly=qy-fk2l*uk2ly
         fk1lz=qz-fk2l*uk2lz

         call vecthph(1,fk1lx,fk1ly,fk1lz,fk1l,th1l(isol),ph1l(isol))

         e1l(isol)=sqrt((hc*fk1l)**2+fm1**2)

        end do

       else if(ivar==3 .or. (ivar>=30 .and. ivar<40)) then

        e1l(1)=varl+fm1
        th1l(1)=thxl
        ph1l(1)=phxl
        e2l(1)=et2l+fm2
        e1l(2)=e1l(1)
        th1l(2)=th1l(1)
        ph1l(2)=ph1l(1)
        e2l(2)=e2l(1)
        ebl=e0l+eal-e1l(1)-e2l(1)
        if(ebl<fmb) return

        fk1l=sqrt((e1l(1)+fm1)*(e1l(1)-fm1))/hc
        fk2l=sqrt((e2l(1)+fm2)*(e2l(1)-fm2))/hc
        fkbl(1)=sqrt((ebl+fmb)*(ebl-fmb))/hc
        fkbl(2)=fkbl(1)

        call vecxyz(1,fk1l,th1l(1),ph1l(1), fk1lx,fk1ly,fk1lz)

        qx=-fk1lx
        qy=-fk1ly
        qz=fk0l+fkal-fk1lz

        call vecthph(1,qx,qy,qz, q,thq,phq)

        cthq=cos(thq*piad)
        sthq=sin(thq*piad)
        cphq=cos(phq*piad)
        sphq=sin(phq*piad)

        cth2p=(fk2l**2+q**2-fkbl(1)**2)/(2.0d0*q*fk2l)

        if(abs(cth2p)>1.0d0+epslim) return

        th2p=acos2(cth2p)/piad
        sth2p=sin(th2p*piad)

        cthbp=(q-fk2l*cth2p)/fkbl(1)
        thbp=acos2(cthbp)/piad
        sthbp=sin(thbp*piad)

        faca=sth2p*cthq
        facb=cth2p*sthq
        facc=sth2p

        if(abs(mod(ph2l-90.0d0,180.0d0))<epslim) then
          d2=(faca**2-facb**2)*facc**2*cphq**2*sphq**2+facc**4*sphq**4
          if(d2<0) return
          fnm=-faca*facb*cphq**2
          dnm=faca**2*cphq**2+facc**2*sphq**2

          do isol=1,2
           if(isol==2 .and. abs(d2)<epslim) return
           fugo=1.0d0
           if(isol==1) fugo=-1.0d0
           cph2p=(fnm+fugo*sqrt(d2))/dnm
           if(abs(cph2p)>1.0d0+epslim) cycle

           do iph=1,2
            ph2p=acos2(cph2p)/piad
            if(abs(mod(ph2p,180.0d0))<epslim .and. iph==2) cycle
            if(iph==2) ph2p=360.0d0-ph2p
            sph2p=sin(ph2p*piad)
            cph2p=cos(ph2p*piad)

            fk2lx=(sth2p*cph2p*cthq*cphq+cth2p*sthq*cphq-sth2p*sph2p*sphq)*fk2l
            fk2ly=(sth2p*cph2p*cthq*sphq+cth2p*sthq*sphq+sth2p*sph2p*cphq)*fk2l
            fk2lz=(-sth2p*cph2p*sthq+cth2p*cthq)*fk2l

            call vecthph(1,fk2lx,fk2ly,fk2lz, fk2lw,th2lw,ph2lw)

            if(abs(fk2l-fk2lw)>epslim) stop 111
            if(abs(sin(ph2lw*piad)-sin(ph2l*piad))>epslim) cycle
            if(abs(cos(ph2lw*piad)-cos(ph2l*piad))>epslim) cycle

            if(iapp(isol)==1 .and. iph==2) stop 112
            iapp(isol)=1
            th2l(isol)=th2lw

            fkblx=qx-fk2lx
            fkbly=qy-fk2ly
            fkblz=qz-fk2lz
            call vecthph(1,fkblx,fkbly,fkblz, fkblw,thbl(isol),phbl(isol))
            if(abs(fkbl(1)-fkblw)>epslim) stop 113

           end do
          end do

         else

          tph2=tan(ph2l*piad)
          d2=(cphq+tph2*sphq)**2 &
     &       *((faca**2-facb**2)*(tph2*cphq-sphq)**2+facc**2*(cphq+tph2*sphq)**2)

          if(d2<0) return
          fnm=-faca*facb*(-tph2*cphq+sphq)**2
          dnm=faca**2*(tph2*cphq-sphq)**2+facc**2*(cphq+tph2*sphq)**2

          do isol=1,2
           if(isol==2 .and. abs(d2)<epslim) return
           fugo=1.0d0
           if(isol==1) fugo=-1.0d0
           cph2p=(fnm+fugo*facc*sqrt(d2))/dnm

           if(abs(cph2p)>1.0d0+epslim) cycle

           do iph=1,2
            ph2p=acos2(cph2p)/piad
            if(abs(mod(ph2p,180.0d0))<epslim .and. iph==2) cycle
            if(iph==2) ph2p=360.0d0-ph2p

            sph2p=sin(ph2p*piad)
            cph2p=cos(ph2p*piad)

            fk2lx=(sth2p*cph2p*cthq*cphq+cth2p*sthq*cphq-sth2p*sph2p*sphq)*fk2l
            fk2ly=(sth2p*cph2p*cthq*sphq+cth2p*sthq*sphq+sth2p*sph2p*cphq)*fk2l
            fk2lz=(-sth2p*cph2p*sthq+cth2p*cthq)*fk2l

            call vecthph(1,fk2lx,fk2ly,fk2lz, fk2lw,th2lw,ph2lw)

            if(abs(fk2l-fk2lw)>epslim) stop 114
            if(abs(sin(ph2lw*piad)-sin(ph2l*piad))>epslim) cycle
            if(abs(cos(ph2lw*piad)-cos(ph2l*piad))>epslim) cycle

            if(iapp(isol)==1 .and. iph==2) stop 115
            iapp(isol)=1
            th2l(isol)=th2lw

            fkblx=qx-fk2lx
            fkbly=qy-fk2ly
            fkblz=qz-fk2lz
            call vecthph(1,fkblx,fkbly,fkblz, fkblw,thbl(isol),phbl(isol))
            if(abs(fkbl(1)-fkblw)>epslim) stop 116

           end do
          end do
        end if

      end if

      nsol=iapp(1)+iapp(2)

      return

end subroutine kinema3bl
!=======================================================================



!=======================================================================
function fkin(fk1c,e1c,th1c,ph1c,fk2c,e2c,th2c,ph2c,fkbc,ebc,thbc,phbc &
     &       ,fktotiz,ivar)

      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      fkin=0.0d0

      fk2cw=fk2c
      e1cw =e1c
      ebcw =ebc

      if(fk2cw<epslim) fk2cw=epslim
      if(e1cw <epslim) e1cw =epslim
      if(ebcw <epslim) ebcw =epslim

      call vecxyz(1,fk2c,th2c,ph2c, fk2cx,fk2cy,fk2cz)

      if(ivar==1 .or. (ivar>9 .and. ivar<20) .or. ivar>39) then
        call vecxyz(1,fk1c,th1c,ph1c, xx,xy,fk1cz)
        xz=fk1cz-fktotiz
        sprox2=xx*fk2cx+xy*fk2cy+xz*fk2cz
        fac1=fk1c*fk2c*e1c*e2c/hc**4
        fac2=e2c/ebcw
        fac3=e2c/ebcw*sprox2/fk2cw**2
        denom=abs(1.0d0+fac2+fac3)
        if(denom<epslim) denom=epslim
        fkin=fac1/denom
       else if(ivar==2 .or. (ivar>19 .and. ivar<30)) then
        call vecxyz(1,fkbc,thbc,phbc, yx,yy,fkbcz)
        yz=fkbcz-fktotiz
        sproy2=yx*fk2cx+yy*fk2cy+yz*fk2cz
        fac1=fk2c*e2c/hc**2
        fac2=e2c/e1cw
        fac3=e2c/e1cw*sproy2/fk2cw**2
        denom=abs(1.0d0+fac2+fac3)
        if(denom<epslim) denom=epslim
        fkin=fac1/denom
       else if(ivar==3 .or. (ivar>29 .and. ivar<40)) then
        call vecxyz(1,fk1c,th1c,ph1c, fk1cx,fk1cy,fk1cz)
        qx=-fk1cx
        qy=-fk1cy
        qz=fktotiz-fk1cz
        denom=abs(-qx*cos(th2c*piad)/sin(th2c*piad)*cos(ph2c*piad) &
     &            -qy*cos(th2c*piad)/sin(th2c*piad)*sin(ph2c*piad)+qz)
        if(denom<epslim) denom=epslim
        fac1=e1c*e2c*ebc*fk1c/hc**6
        fkin=fac1/denom
       else
        stop 1
      end if

      return

end function fkin
!=======================================================================



!=======================================================================
function fjac(e1,e2,eb,e1c,e2c,ebc)
       use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      fjac=0.0d0

      e1cw=e1c
      e2cw=e2c
      ebcw=ebc

      if(e1cw<epslim) e1cw=epslim
      if(e2cw<epslim) e2cw=epslim
      if(ebcw<epslim) ebcw=epslim

      fjac=e1*e2*eb/(e1cw*e2cw*ebcw)

      return

end function fjac
!=======================================================================



!=======================================================================
function acos2(costh)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      if(abs(costh)>1.0d0+epslim) then
        write(*,*) 'ERROR: costh has a strange value.'
        write(*,'(a,f0.4)') ' costh=',costh
        stop 1
      end if
      if(costh> 1.0d0) costh= 1.0d0
      if(costh<-1.0d0) costh=-1.0d0
      acos2=acos(costh)

      return

end function acos2
!=======================================================================



!=======================================================================
subroutine vecthph(ictr,vecx,vecy,vecz, vec,theta,phi)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----ictr=1: theta and phi are given in the unit of degree
!          2: theta and phi are given in the unit of radian

!-----------------------------------------------------------------------

      vec=sqrt(vecx**2+vecy**2+vecz**2)

      if(vec<epslim) then
        vec=0.0d0
        theta=0.0d0
        phi=0.0d0
        return
      end if

      costh=vecz/vec
      theta=acos2(costh)

      vecb=sqrt(vecx**2+vecy**2)

      if(vecb<epslim) then
        cosph=1.0d0
        phi=0.0d0
       else
        cosph=vecx/vecb
        phi=acos2(cosph)
        if(vecy<-epslim) phi=2.0d0*pi-phi
      end if

      if(ictr==1) then
        theta=theta/piad
        phi=phi/piad
      end if

      return

end subroutine vecthph
!=======================================================================



!=======================================================================
subroutine vecxyz(ictr,vec,theta,phi,vecx,vecy,vecz)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----ictr=1: theta and phi are given in the unit of degree
!          2: theta and phi are given in the unit of radian

!-----------------------------------------------------------------------

      if(ictr==1) then
        tht=theta*piad
        ph=phi*piad
       else if(ictr==2) then
        tht=theta
        ph=phi
       else
        write(*,*) 'ERROR: ictr must be 1 or 2'
        write(*,'(a,i0)') 'ictr=',ictr
        stop 1
      end if

      vecx=vec*sin(tht)*cos(ph)
      vecy=vec*sin(tht)*sin(ph)
      vecz=vec*cos(tht)

      return

end subroutine vecxyz
!=======================================================================



!=======================================================================
subroutine ltrgen(fkx,fky,fkz,e,betx,bety,betz,bet,gam,fkpx,fkpy,fkpz,ep)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----lorentz transformation for beta in general direction

!-----------------------------------------------------------------------

      fac0=(gam-1.0d0)/bet**2
      prod=fkx*betx+fky*bety+fkz*betz
      fac1=fac0*prod
      fac2=gam/hc*e
      fkpx=fkx+fac1*betx-fac2*betx
      fkpy=fky+fac1*bety-fac2*bety
      fkpz=fkz+fac1*betz-fac2*betz
      ep=gam*(e-prod*hc)

      return

end subroutine ltrgen
!=======================================================================



!=======================================================================
subroutine ltrz(ictr1,var,ictr2,tht,phi,bet,gam,fm,ep,fkp,thp,php)
      use consts
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----lorentz transformation for beta // z

!-----ictr1=1: var is the enegy
!           2: var is the momentum (in the unit of hbar)

!-----ictr2=1: theta and phi are given in the unit of degree
!           2: theta and phi are given in the unit of radian

!-----------------------------------------------------------------------

      if(ictr1==1) then
        e=var
        fk=sqrt((e+fm)*(e-fm))/hc
       else if(ictr1==2) then
        fk=var
        e=sqrt((hc*fk)**2+fm**2)
       else
        write(*,*) 'ERROR: ictr1 must be 1 or 2'
        write(*,'(a,i0)') 'ictr1=',ictr1
        stop 1
      end if

      call vecxyz(ictr2,fk,tht,phi,fkx,fky,fkz)

      fkpz=gam*(fkz-bet*e/hc)
      ep=gam*(e-bet*hc*fkz)

      call vecthph(ictr2,fkx,fky,fkpz, fkp,thp,php)
      php=phi

      return

end subroutine ltrz
!=======================================================================



!=======================================================================
subroutine input(limfs,ions,ifrm,imir,ical                   &
     &          ,zp,ap,za,aa,ikin,elab,ictrein               &
     &          ,ish,ebind,zsp,asp,betasp,ictrm              &
     &          ,fj,ls,sfac,nod,ibmc,rc,ictrc,a0c,rcl,ictrcl &
     &          ,ibms,v0ls,rs,ictrs,as,lmax                  &
     &          ,ivar,iex,fkncut,ixunt,kunt                  &
     &          ,ivvarl,varlmin,varlmax,dvarl                &
     &          ,ivthxl,thxlmin,thxlmax,dthxl                &
     &          ,ivphxl,phxlmin,phxlmax,dphxl                &
     &          ,ivet2l,et2lmin,et2lmax,det2l                &
     &          ,ivph2l,ph2lmin,ph2lmax,dph2l                &
     &          ,ielm,kinelm,ielmedg                         &
     &          ,rmax,dr,ngth,ngph,ngk1,ngph1,nrgmax         &
     &          ,ipot,facv,facw,facvs,facws,beta,ims,iedg    &
     &          ,detbl,itmdcal,icpot,ilscal)
      use consts
      use dims
      use kibmod
      use nntbl,only:ionsh
      use array,only:lauto
      implicit real*8 (a-h,o-z)
!=======================================================================
      character*2 namep,namet,namesp,nameb
      character*50 comment
      character*5 strlmax(3)

      integer ipot(3)
      real*8 facv(3),facw(3),facvs(3),facws(3),beta(3)
      integer ims(3),iedg(3),lmax(3)
      integer icpot(3)

!-----------------------------------------------------------------------

      kiban=5

!-----default values
      limfs=1000
      ions=0
      ifrm=0
      imir=0
      ical=1
      ikin=0
      ictrein=0
      betasp=0.0d0
      ictrm=1
      iex=0
      fkncut=2.0d0
      ixunt=0
      kunt=1
      kibtmd=0
      kiblg=0
      kibpx=0
      kibtr=0
      kibbs=0
      kibtl=0
      ionsh=1
      kinelm=0
      ielmedg=1
      ngb=30
      ngth=40
      ngph=40
      ngk1=30
      ngph1=1
      facv(:)=1.0d0
      facw(:)=1.0d0
      facvs(:)=1.0d0
      facws(:)=1.0d0
      beta(:)=0.0d0
      ims(:)=0
      iedg(:)=1

      read(kiban,500) comment
  500 format(a50)
      read(kiban,501) limfs,ions,ifrm,imir,ical
  501 format(5i5)
      if(limfs==0) limfs=1000000
      read(kiban,502) zp,ap,za,aa
  502 format(f5.0,f10.0,f5.0,f10.0)
      read(kiban,506) ikin,elab,ictrein
  503 format(i5,4f10.0)
      read(kiban,504) ish,ebind,zsp,asp,betasp,ictrm
  504 format(i5,f10.0,f5.0,2f10.0,i5)
      read(kiban,505) fj,fl,sfac,nod,kibbs
  505 format(2f5.0,f10.0,2i5)
      read(kiban,506) ibmc,rc,ictrc,a0c,rcl,ictrcl
  506 format(i5,f10.0,i5,2f10.0,i5)
      read(kiban,507) ibms,v0ls,rs,ictrs,as
  507 format(i5,2f10.0,i5,f10.0)
      read(kiban,513) lmax(1),lmax(2),lmax(3)
      read(kiban,508) ivar,iex,fkncut,ixunt,kunt
  508 format(2i5,f5.0,2i5)
  511 format(7i5)
  512 format(5i5)
  513 format(3i5)
      read(kiban,503) ivvarl,varlmin,varlmax,dvarl
      read(kiban,503) ivthxl,thxlmin,thxlmax,dthxl
      read(kiban,503) ivphxl,phxlmin,phxlmax,dphxl
      read(kiban,503) ivet2l,et2lmin,et2lmax,det2l
      read(kiban,503) ivph2l,ph2lmin,ph2lmax,dph2l
      read(kiban,511) kibtbl,kibout,kibtmd,kiblg,kibpx,kibtr,kibtl
      read(kiban,512) ielm,kibelm,ionsh,kinelm,ielmedg
      read(kiban,509) rmax,dr,ngb,ngth,ngph,ngk1,ngph1
  509 format(2f10.0,5i5)
      read(kiban,510) ipot(1),facv(1),facw(1),facvs(1),facws(1),beta(1),ims(1),iedg(1)
      read(kiban,510) ipot(2),facv(2),facw(2),facvs(2),facws(2),beta(2),ims(2),iedg(2)
      read(kiban,510) ipot(3),facv(3),facw(3),facvs(3),facws(3),beta(3),ims(3),iedg(3)
  510 format(i5,4f5.0,f10.0,2i5)
      detbl=5.d0


!-------------------------------
!---  input condition check  ---
!-------------------------------

       if(ifrm<0 .or. ifrm>2) then
         write(*,*) 'ERROR: ifrm must be 0--2'
         write(*,'(a,i0)') ' ifrm=',ifrm
         stop 1
       end if

       if(zp<0.0d0) then
         write(*,*) 'ERROR: zp cannot be negative'
         write(*,'(a,f0.4)') ' zp=',zp
         stop 1
       end if

       if(ap<=0.0d0) then
         write(*,*) 'ERROR: ap must be positive'
         write(*,'(a,f0.4)') ' ap=',ap
         stop 1
       end if

       if(ap<zp) then
         write(*,*) 'ERROR: ap smaller than zp'
         write(*,'(a,f0.4,a,f0.4)') ' zp=',zp,', ap=',ap
         stop 1
       end if

       if(za<0.0d0) then
         write(*,*) 'ERROR: za cannot be negative'
         write(*,'(a,f0.4)') ' za=',za
         stop 1
       end if

       if(aa<=0.0d0) then
         write(*,*) 'ERROR: aa must be positive'
         write(*,'(a,f0.4)') ' aa=',aa
         stop 1
       end if

       if(aa<za) then
         write(*,*) 'ERROR: aa smaller than za'
         write(*,'(a,f0.4,a,f0.4)') ' za=',za,', aa=',aa
         stop 1
       end if

       if(zsp<0.0d0) then
         write(*,*) 'ERROR: zsp cannot be negative'
         write(*,'(a,f0.4)') ' zsp=',zsp
         stop 1
       end if

       if(asp<=0.0d0) then
         write(*,*) 'ERROR: asp must be positive'
         write(*,'(a,f0.4)') ' asp=',asp
         stop 1
       end if

       if(asp<zsp) then
         write(*,*) 'ERROR: asp smaller than zsp'
         write(*,'(a,f0.4,a,f0.4)') ' zsp=',zsp,', asp=',asp
         stop 1
       end if

       if(za-zsp<0.0d0) then
         write(*,*) 'ERROR: za-zsp cannot be negative'
         write(*,'(a,f0.4)') ' za-zsp=',za-zsp
         stop 1
       end if

       if(aa-asp<=0.0d0) then
         write(*,*) 'ERROR: aa-asp must be positive'
         write(*,'(a,f0.4)') ' aa-asp=',aa-asp
         stop 1
       end if

       if(ish<0 .or. (ish>1 .and. ish<10)) then
         write(*,*) 'ERROR: ish must be 0, 1 or larger than 9'
         write(*,'(a,i0)') ' ish=',ish
         stop 1
       end if

       if(ictrm<1 .or. ictrm>2) then
         write(*,*) 'ERROR: ictrm must be 1 or 2'
         write(*,'(a,i0)') ' ictrm=',ictrm
         stop 1
       end if

       if(ibmc/=1) then
         if(rc<=0.0d0) then
           write(*,*) 'ERROR: rc must be positive'
           write(*,'(a,f0.4)') ' rc=',rc
           stop 1
         end if

         if(ictrc<0 .or. ictrc>3) then
           write(*,*) 'ERROR: ictrc must be 0--3'
           write(*,'(a,i0)') ' ictrc=',ictrc
           stop 1
         end if

         if(a0c<=0.0d0) then
           write(*,*) 'ERROR: a0c must be positive'
           write(*,'(a,f0.4)') ' a0c=',a0c
           stop 1
         end if

         if(rcl<=0.0d0) then
           write(*,*) 'ERROR: rcl must be positive'
           write(*,'(a,f0.4)') ' rcl=',rcl
           stop 1
         end if

         if(ictrcl<0 .or. ictrcl>3) then
           write(*,*) 'ERROR: ictrcl must be 0--3'
           write(*,'(a,i0)') ' ictrcl=',ictrcl
           stop 1
         end if
       end if

       if(ibms/=1) then
         if(rs<=0.0d0) then
           write(*,*) 'ERROR: rs must be positive'
           write(*,'(a,f0.4)') ' rs=',rs
           stop 1
         end if

         if(ictrs<0 .or. ictrs>3) then
           write(*,*) 'ERROR: ictrs must be 0--3'
           write(*,'(a,i0)') ' ictrs=',ictrs
           stop 1
         end if

         if(as<=0.0d0) then
           write(*,*) 'ERROR: as must be positive'
           write(*,'(a,f0.4)') ' as=',as
           stop 1
         end if
       end if

       if((ibmc==1 .or. ibms==1) .and. nint(asp)/=1) then
         write(*,*) 'ERROR: ibmc/ibms=1 is not allowed for non-nucleon knockout'
         write(*,'(a,i0,a,i0,a,f0.4)') ' ibmc=',ibmc,' ibms=',ibms,' asp=',asp
         stop 1
       end if

       lov=0
       do i=1,3
        if(lmax(i)>=0 .and. lmax(i)+1>ldim) then
          lov=1
          write(*,*) 'ERROR: lmax for particle i > ldim'
          write(*,'(a,i0,a,i0,a,i0)') ' lmax+1=',lmax(i)+1,', i=',i-1,', ldim=',ldim
        end if
       end do
       if(lov==1) stop 1

       lov=0
       do i=1,3
        if(lmax(i)<0 .and. abs(lmax(i))+1>ldim) then
          lov=1
          write(*,*) 'ERROR: |lmax| for particle i > ldim'
          write(*,'(a,i0,a,i0,a,i0)') &
     &      ' |lmax|+1=',abs(lmax(i))+1,', i=',i-1,', ldim=',ldim
        end if
       end do
       if(lov==1) stop 1

       if(ivar<=0 .or. (ivar>3 .and. ivar<9)) then
         write(*,*) 'ERROR: ivar must be 1, 2, 3, 9, or >9'
         write(*,'(a,i0)') ' ivar=',ivar
         stop 1
       end if

      if(kunt<0 .or. kunt>2) then
        write(*,*) 'ERROR: kunt must be 0--2'
        write(*,'(a,i0)') ' kunt=',kunt
        stop 1
      end if

      if(ivvarl/=0) then
        if(varlmin>varlmax) then
          write(*,*) 'ERROR: varlmin > varlmax'
          write(*,'(a,f0.4,a,f0.4)') ' varlmin=',varlmin,', varlmax=',varlmax
          stop 1
        end if
        if(dvarl<=0.0d0) then
          write(*,*) 'ERROR: dvarl must be positive'
          write(*,'(a,f0.4)') ' dvarl=',dvarl
          stop 1
        end if
      end if

      if(ivthxl/=0) then
        if(thxlmin>thxlmax) then
          write(*,*) 'ERROR: thxlmin > thxlmax'
          write(*,'(a,f0.4,a,f0.4)') ' thxlmin=',thxlmin,', thxlmax=',thxlmax
          stop 1
        end if
        if(dthxl<=0.0d0) then
          write(*,*) 'ERROR: dthxl must be positive'
          write(*,'(a,f0.4)') ' dthxl=',dthxl
          stop 1
        end if
      end if

      if(ivphxl/=0) then
        if(phxlmin>phxlmax) then
          write(*,*) 'ERROR: phxlmin > phxlmax'
          write(*,'(a,f0.4,a,f0.4)') ' phxlmin=',phxlmin,', phxlmax=',phxlmax
          stop 1
        end if
        if(dphxl<=0.0d0) then
          write(*,*) 'ERROR: dphxl must be positive'
          write(*,'(a,f0.4)') ' dphxl=',dphxl
          stop 1
        end if
      end if

      if(ivet2l/=0) then
        if(et2lmin>et2lmax) then
          write(*,*) 'ERROR: et2lmin > et2lmax'
          write(*,'(a,f0.4,a,f0.4)') ' et2lmin=',et2lmin,', et2lmax=',et2lmax
          stop 1
        end if
        if(det2l<=0.0d0) then
          write(*,*) 'ERROR: det2l must be positive'
          write(*,'(a,f0.4)') ' det2l=',det2l
          stop 1
        end if
      end if

      if(ivph2l/=0) then
        if(ph2lmin>ph2lmax) then
          write(*,*) 'ERROR: ph2lmin > ph2lmax'
          write(*,'(a,f0.4,a,f0.4)') ' ph2lmin=',ph2lmin,', ph2lmax=',ph2lmax
          stop 1
        end if
        if(dph2l<=0.0d0) then
          write(*,*) 'ERROR: dph2l must be positive'
          write(*,'(a,f0.4)') ' dph2l=',dph2l
          stop 1
        end if
      end if

      if(kibtbl<=0) then
        write(*,*) 'ERROR: kibtbl must be an positive integer'
        write(*,'(a,i0)') ' kibtbl=',kibtbl
        stop 1
      end if

      if(kibout<=0) then
        write(*,*) 'ERROR: kibout must be an positive integer'
        write(*,'(a,i0)') ' kibout=',kibout
        stop 1
      end if

      if(ivar==9) then
        ierrmd=0
        nkbzmax=nint((varlmax-varlmin)/dvarl)+1
        nkbbmax=nint((thxlmax-thxlmin)/dthxl)+1
        if(ivvarl==0) nkbzmax=1
        if(ivthxl==0) nkbbmax=1
        if(kibpx>0 .and. nkbbmax<4) then
          write(*,*) 'ERROR: more than 3 points of K_Bb needed to calculate ds/dK_Bx'
          write(*,'(a,i0)') ' # of K_Bb to be calculated=',nkbbmax
          ierrmd=1
        end if
        if(kibtl>0 .and. (nkbbmax<4 .or. nkbzmax<4)) then
          write(*,*) 'ERROR: more than 3 points of K_Bz and K_Bb needed to calculate' &
     &              ,' ds/dK_Btot'
          write(*,'(a,i0)') ' # of K_Bz to be calculated=',nkbzmax
          write(*,'(a,i0)') ' # of K_Bb to be calculated=',nkbbmax
          ierrmd=1
        end if
        if(ielm==4 .and. ngph1/=1) then 
          write(*,'(a)') 'ngph1 must be 1 when IVAR=9 and IELM=4'
          ierrmd=1
        end if
        if(ierrmd==1) stop 1
      end if

      if(ielm<0 .or. ielm==1 .or. ielm==2 .or.ielm>4) then
        write(*,*) 'ERROR: ielm must be 0, 3, or 4'
        write(*,'(a,i0)') ' ielm=',ielm
        stop 1
      end if

      if(ionsh<1 .or. ionsh>4) then
        write(*,*) 'ERROR: ionsh must be 1--4'
        write(*,'(a,i0)') ' ionsh=',ionsh
        stop 1
      end if

       if(rmax<=0.0d0) then
         write(*,*) 'ERROR: rmax must be positive'
         write(*,'(a,f0.4)') ' rmax=',rmax
         stop 1
       end if

       if(dr<=0.0d0) then
         write(*,*) 'ERROR: dr must be positive'
         write(*,'(a,f0.4)') ' dr=',dr
         stop 1
       end if

      if(max(ngb,ngth,ngph,ngk1,ngph1)>ngdim) then
        write(*,*) 'ERROR: one of ng-x > ngdim'
        write(*,'(a,i0,a,i0)') ' ng-x=',max(ngb,ngth,ngph,ngk1,ngph1),', ngdim=',ngdim
        stop 1
      end if

      if(ngb*ngth*ngph==0) then
        write(*,*) 'ERROR: one of ng-b,th,ph is 0'
        write(*,'(a,3i3)') ' ngb ngth ngph=',ngb,ngth,ngph
        stop 1
      end if

      if(ivar==9 .and. ngk1*ngph1==0) then
        write(*,*) 'ERROR: one of ng-k1,ph1 is 0 with ivar=9'
        write(*,'(a,2i3)') ' ngk1 ngph1=',ngk1,ngph1
        stop 1
      end if

      if(betasp<-epslim .and. ish<10) then
        write(*,*) 'ERROR: betasp<0 is allowed only when ish>9'
        write(*,'(a,f0.4,a,i0)') 'betasp=',betasp,', ish=',ish
        stop 1
      end if

      if(ish>9 .and. betasp>epslim .and. ibmc/=1) then
        write(*,*) 'ERROR: ibmc must be 1 when ish>9 and betasp>0'
        write(*,'(a,i0,a,i0,a,f0.4)') 'ibmc=',ibmc,', ish=',ish,', betasp=',betasp
        stop 1
      end if

      if(beta(1)<-epslim .and. ipot(1)<10 .and. ipot(1)/=0) then
        write(*,*) 'ERROR: beta(1)<0 is allowed only when ipot(1)>9 for particle 0'
        write(*,'(a,f0.4,a,i0)') 'beta(1)=',beta(1),', ipot(1)=',ipot(1)
        stop 1
      end if

      if(beta(2)<-epslim .and. ipot(2)<10 .and. ipot(2)/=0) then
        write(*,*) 'ERROR: beta(2)<0 is allowed only when ipot(2)>9 for particle 1'
        write(*,'(a,f0.4,a,i0)') 'beta(2)=',beta(2),', ipot(2)=',ipot(2)
        stop 1
      end if

      if(beta(3)<-epslim .and. ipot(3)<10 .and. ipot(3)/=0) then
        write(*,*) 'ERROR: beta(3)<0 is allowed only when ipot(3)>9 for particle 2'
        write(*,'(a,f0.4,a,i0)') 'beta(3)=',beta(3),', ipot(3)=',ipot(3)
        stop 1
      end if

      ivsum=ivvarl+ivthxl+ivphxl+ivet2l+ivph2l
      if(ivsum>1 .and. kibtmd>=1) then
        write(*,*) 'ERROR: tmd output is allowed only when single or no degree of' &
     &             ,' freedom in kinematics is varied'
        stop 1
      end if
      if(ielm==4 .and. kibtmd>=1) then
        write(*,*) 'ERROR: tmd output is not supported for IELM = 4'
        stop 1
      end if

!-----restriction of the current version----------------------
      if(ap>1.1d0) then
         write(*,*) 'ERROR: the probe must be p or n'
         stop 1
       end if

      if(iex==1 .and. (ielm==4 .or. ielm==6)) then
        write(*,*) 'ERROR: iex=1 is not supported when IELM=4 or 6'
        stop 1
      end if

      if(asp>1.1d0 .and. ielm/=3) then
         write(*,*) 'ERROR: NN scattering assumed for ielm/=3'
         stop 1
       end if

      if(asp>1.1d0 .and. ipot(3)==1) then
         write(*,*) 'ERROR: KD potential not applicable to particle not a nucleon'
         stop 1
       end if

      if(elab>200.0d0 .and. ipot(1)==1) then
         write(*,*) 'ERROR: KD potential can be used for energy below 200 MeV'
         stop 1
       end if

      if(ipot(1)/=0 .and. ipot(1)/=1 .and. ipot(1)<10) then
         write(*,*) 'ERROR: ipot must be 0, 1 or larger than 9 for particle 0'
         write(*,'(a,i0)') ' ipot(1)=',ipot(1)
         stop 1
       end if

      if(ipot(2)/=0 .and. ipot(2)/=1 .and. ipot(2)<10) then
         write(*,*) 'ERROR: ipot must be 0, 1 or larger than 9 for particle 1'
         write(*,'(a,i0)') ' ipot(2)=',ipot(2)
         stop 1
       end if

      if(ipot(3)/=0 .and. ipot(3)/=1 .and. ipot(3)<10) then
         write(*,*) 'ERROR: ipot must be 0, 1 or larger than 9 for particle 2'
         write(*,'(a,i0)') ' ipot(3)=',ipot(3)
         stop 1
       end if

      if(ielm==4) then
        if(ivphxl+ivph2l>0) then
          write(*,*) 'ERROR: phi cannot be varied when IELM=4'
          stop 1
        end if
        if(abs(phxlmin)>1.d-3.and.abs(phxlmin-180.d0)>1.d-3) then
          write(*,*) 'ERROR: kinamatics must be coplanar when IELM=4'
          stop 1
        end if
        if(abs(ph2lmin)>1.d-3.and.abs(ph2lmin-180.d0)>1.d-3) then
          write(*,*) 'ERROR: kinamatics must be coplanar when IELM=4'
          stop 1
        end if
        if(abs(phxlmin-ph2lmin)-180.d0>1.d-3) then
          write(*,*) 'ERROR: kinamatics must be coplanar when IELM=4'
          stop 1
        end if
        if(ap > 1.1d0 .or. asp > 1.1d0) then
          write(*,*) 'ERROR: only (p,pN) reaction is allowed when IELM=4 or 6'
          stop 1
        end if
      end if
!-------------------------------------------------------------

      ichng1=0

      ichng2=0
      if((ical==0 .or. ivar==9) .and. kinelm==1) then
        ichng2=1
        kinelmw=kinelm
        kinelm=0
      end if

      ichng3=0
      if(ical==0 .and. ivar==9) then
        ichng3=1
        icalw=ical
        ical=1
      end if

      ichng4=0
      if(ivar==9) then
        if(ikin==0 .and. ifrm/=0) then
          ichng4=1
          ifrmw=ifrm
          ifrm=0
         else if(ikin==1 .and. ifrm/=2) then
          ichng4=1
          ifrmw=ifrm
          ifrm=2
        end if
      end if


!----------------------------------------
!---  names of projectile and target  ---
!----------------------------------------

      call elemnt(ap,zp,massp,namep)
      call elemnt(aa,za,masst,namet)
      call elemnt(asp,zsp,masssp,namesp)
      ab=aa-asp
      zb=za-zsp
      call elemnt(ab,zb,massb,nameb)


!-------------------------------
!---  output for the inputs  ---
!-------------------------------

      write(kibout,600)
  600 format(5x,'*****************************************************************' &
     &      /5x,'*****                                                       *****' &
     &      /5x,'*****                 program pikoe ver.1.0                 *****' &
     &      /5x,'*****                                                       *****' &
     &      /5x,'*****************************************************************')

      if(iex==0) then
        write(kibout,601) massp,namep,masst,namet,nint(elab) &
     &                   ,massp,namep,masssp,namesp,massb,nameb
       else if(iex==1) then
        write(kibout,601) massp,namep,masst,namet,nint(elab) &
     &                   ,masssp,namesp,massp,namep,massb,nameb
      end if
  601 format(//3x,'============================================================='   &
     &        /8x,'reaction:',3x,i3,a2,' + ',i3,a2,5x,'at ',i5,'MeV/nucleon'        &
     &       //8x,'particle 1:',i3,a2,2x,'particle 2:',i3,a2,2x,'particle B:',i3,a2 &
     &        /3x,'=============================================================')

      write(kibout,602) comment
  602 format(/5x,'user''s comment:'/10x,a50)

      write(kibout,603) limfs,ions,ifrm,imir,ical, zp,ap,za,aa &
     &                 ,ikin,elab,ictrein
  603 format(/3x,'-- general inputs --'                              &
     &       /8x,'limfs ions ifrm imir ical        :',5i5            &
     &       /8x,'zp ap za aa                      :',2(f10.2,f10.5) &
     &       /8x,'ikin elab ictrein                :',i5,f10.5,i5)

      if(ish==0) then
        write(kibout,604) ish      ,zsp,asp,betasp,ictrm,fj,fl,sfac,nod,kibbs
       else
        write(kibout,605) ish,ebind,zsp,asp,betasp,ictrm,fj,fl,sfac,nod,kibbs
      end if
  604 format( 8x,'ish ebind zsp asp betasp ictrm   :',i5,' not-fixed' &
     &                                               ,f10.2,2f10.5,i5 &
     &       /8x,'fj  fl sfac  nod kibbs           :',3f10.2,2i5)
  605 format( 8x,'ish ebind zsp asp betasp ictrm   :',i5,f10.5        &
     &                                               ,f10.2,2f10.5,i5 &
     &       /8x,'fj  fl sfac  nod kibbs           :',3f10.2,2i5)

      if(ibmc==1) then
        if(ish==0) then
          write(kibout,606) 1.27,1,0.67,1.27,1
         else
          write(kibout,607) 1.27,1,0.67,1.27,1
        end if
       else
        if(ish==0) then
          write(kibout,608) ebind,rc,ictrc,a0c,rcl,ictrcl
         else
          write(kibout,609)       rc,ictrc,a0c,rcl,ictrcl
        end if
      end if
  606 format(/8x,'-vce- Bohr-Mottelson'                         &
     &       /8x,'    v0ce rc ictrc a0c rcl ictrcl : default  ' &
     &                                              ,f10.5,i5,2f10.5,i5)
  607 format(/8x,'-vce- Bohr-Mottelson'                         &
     &       /8x,'    v0ce rc ictrc a0c rcl ictrcl : not-fixed' &
     &                                              ,f10.5,i5,2f10.5,i5)
  608 format(/8x,'-vce-'                                        &
     &       /8x,'    v0ce rc ictrc a0c rcl ictrcl :',2f10.5,i5 &
     &                                              ,2f10.5,i5)
  609 format(/8x,'-vce-'                                        &
     &       /8x,'    v0ce rc ictrc a0c rcl ictrcl : not-fixed' &
     &                                              ,f10.5,i5,2f10.5,i5)

      if(ibms==1) then
        write(kibout,610) 1.27,1,0.67
       else
        write(kibout,611) v0ls,rs,ictrs,as
      end if
  610 format( 8x,'-vls- Bohr-Mottelson'                         &
     &       /8x,'    v0ls rs ictrs as             : default  ' &
     &                                              ,f10.5,i5,f10.5)
  611 format( 8x,'-vls-' &
     &       /8x,'    v0ls rs ictrs as             :',2f10.5,i5,f10.5)

      strlmax(:)=' auto'
      do i=1,3
       if(lmax(i)>=0) write(strlmax(i),'(i5)') lmax(i)
      end do
      write(kibout,612) strlmax(1),strlmax(2),strlmax(3) &
     &                 ,ivar,iex,fkncut,ixunt,kunt
  612 format(/8x,'lmax0 lmax1 lmax2                :',3a5 &
     &       /8x,'ivar iex fkncut ixunt kunt       :',2i5,f10.5,2i5)

      if(ivar<10) then
        if(ivar==1) then
          write(kibout,613) ivvarl,varlmin,varlmax,dvarl &
     &                     ,ivthxl,thxlmin,thxlmax,dthxl &
     &                     ,ivphxl,phxlmin,phxlmax,dphxl
  613     format(/8x,'ivt1l  t1lmin  t1lmax  dt1l      :',i5,3f10.5 &
     &           /8x,'ivth1l th1lmin th1lmax dth1l     :',i5,3f10.5 &
     &           /8x,'ivph1l ph1lmin ph1lmax dph1l     :',i5,3f10.5)
         else if(ivar==2) then
          write(kibout,614) ivvarl,varlmin,varlmax,dvarl &
     &                     ,ivthxl,thxlmin,thxlmax,dthxl &
     &                     ,ivphxl,phxlmin,phxlmax,dphxl
  614     format(/8x,'ivkbl  kblmin  kblmax  dkbl      :',i5,3f10.5 &
     &           /8x,'ivthbl thblmin thblmax dthbl     :',i5,3f10.5 &
     &           /8x,'ivphbl phblmin phblmax dphbl     :',i5,3f10.5)
         else if(ivar==9) then
          write(kibout,615) ivvarl,varlmin,varlmax,dvarl &
     &                     ,ivthxl,thxlmin,thxlmax,dthxl
  615     format(/8x,'ivkbaz kbazmin kbazmax dkba      :',i5,3f10.5 &
     &           /8x,'ivkbab kbabmin kbabnax dkba      :',i5,3f10.5)
        end if
        if(ivar<=2)                                      &
     &    write(kibout,616) ivet2l,et2lmin,et2lmax,det2l &
     &                     ,ivph2l,ph2lmin,ph2lmax,dph2l
  616   format( 8x,'ivth2l th2lmin th2lmax dth2l     :',i5,3f10.5 &
     &         /8x,'ivph2l ph2lmin ph2lmax dph2l     :',i5,3f10.5)

        if(ivar==3)                                      &
     &    write(kibout,619) ivet2l,et2lmin,et2lmax,det2l &
     &                     ,ivph2l,ph2lmin,ph2lmax,dph2l
  619   format( 8x,'ivt2l  t2lmin  t2lmax  dt2l      :',i5,3f10.5 &
     &         /8x,'ivph2l ph2lmin ph2lmax dph2l     :',i5,3f10.5)

       else
        write(kibout,617)
  617   format(/8x,'kinematical condition is given in external file')
      end if
      write(kibout,618) kibtbl,kibout,kibtmd,kiblg,kibpx,kibtr,kibtl &
     &            ,ielm,kibelm,ionsh,kinelm,ielmedg                  &
     &            ,rmax,dr,ngb,ngth,ngph,ngk1,ngph1                  &
     &            ,ipot(1),facv(1),facw(1),facvs(1),facws(1)         &
     &            ,beta(1),ims(1),iedg(1)                            &
     &            ,ipot(2),facv(2),facw(2),facvs(2),facws(2)         &
     &            ,beta(2),ims(2),iedg(2)                            &
     &            ,ipot(3),facv(3),facw(3),facvs(3),facws(3)         &
     &            ,beta(3),ims(3),iedg(3)
  618 format(/8x,'kib: tbl out tmd  lg  px  tr  tl :',7i5          &
     &       /8x,'ielm kibelm ionsh kinelm ielmedg :',5i5          &
     &       /8x,'rmax  dr                         :',2f10.5       &
     &       /8x,'ng   for  r(b)  th  ph  k1  ph1  :',5i5          &
     &      //8x,'ipot0 facv0 facw0 facvs0 facws0  :',i5,4f10.5    &
     &       /8x,'      beta0 ims0  iedg0          :',5x,f10.5,2i5 &
     &       /8x,'ipot1 facv1 facw1 facvs1 facws1  :',i5,4f10.5    &
     &       /8x,'      beta1 ims1  iedg1          :',5x,f10.5,2i5 &
     &       /8x,'ipot2 facv2 facw2 facvs2 facws2  :',i5,4f10.5    &
     &       /8x,'      beta2 ims2  iedg2          :',5x,f10.5,2i5,/)

      if(ichng1+ichng2+ichng3+ichng4/=0) write(kibout,651)
  651 format(/3x,'-- information on the change(s) of the inputs --')
      if(ichng1==1) write(kibout,652) ionshw,ionsh
  652 format(8x,'ionsh   :'i5,'  -->',i5)
      if(ichng2==1) write(kibout,653) kinelmw,kinelm
  653 format(8x,'kinelm  :'i5,'  -->',i5)
      if(ichng3==1) write(kibout,654) icalw,ical
  654 format(8x,'ical    :'i5,'  -->',i5)
      if(ichng4==1) write(kibout,655) ifrmw,ifrm
  655 format(8x,'ifrm    :'i5,'  -->',i5)

      if(ivvarl==0) varlmax=varlmin
      if(ivthxl==0) thxlmax=thxlmin
      if(ivphxl==0) phxlmax=phxlmin
      if(ivet2l==0) et2lmax=et2lmin
      if(ivph2l==0) ph2lmax=ph2lmin

      nrgmax=ngb

!----flag for option
      itmdcal=0
      if(ivar <= 3) then
        if(kibtmd>=1) itmdcal=1
      end if

      ichngx=0
      nocou=0
      do i=1,3
       if(facv(i)<0.d0) then
         nocou=1
         icpot(i)=0
         facv(i)=abs(facv(i))
         write(kibout,'(/3x,a)') '-- changes from the input (distortion/spin) --'
         ichngx=1
         write(kibout,'(8x,a,i1,a,i1,a)') &
     &     'facv',i-1,' is negetive: Coulomb potential of particle ',i-1,' is turned off'
        else
         icpot(i)=1
       end if
      end do
      if(nocou==1) write(kibout,*)

      if(ielm/=4 .and. ielm/=6) then
        facvs(:)=0.d0
        facws(:)=0.d0
        if(ichngx==0) then
          write(kibout,'(/3x,a)') '-- changes from the input (distortion/spin) --'
          ichngx=1
        end if
        write(kibout,'(8x,a,i1,a)') &
     &    'IELM=',ielm,': Spin-orbit potential of distorted waves are neglected'
        write(kibout,'(12x,a,/)') '-> All facvs = facws = 0.'
      end if

      faclssum=sum(facvs(:))+sum(facws(:))
      if(abs(faclssum)<epslim) then
        ilscal=0
        if(ichngx==0) then
          write(kibout,'(/3x,a)') '-- changes from the input (distortion/spin) --'
          ichngx=1
        end if
        write(kibout,'(8x,a)') &
     &    'All facvs and facws are zero'
        write(kibout,'(12x,a,/)') '-> The spin degree of freedom is neglected.'
       else
        ilscal=1
      end if

      do i=1,3
       if(lmax(i) < 0) then
         lmax(i)=abs(lmax(i))
         lauto(i)=1
        else
         lauto(i)=0
       end if
      end do
!----flag for option end

      ls=nint(fl)

      return

end subroutine input
!=======================================================================



!=======================================================================
subroutine wspotgen(aa,v0,r0r,ar,w0,r0i,ai,wd,r0d,ad &
     &             ,vso,r0rs,ars,wso,r0is,ais        &
     &             ,potv,potw,potd,potrs,potis,nrmax,dr)
      use consts
      use dims
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 potv(nrdim),potw(nrdim),potd(nrdim),potrs(nrdim),potis(nrdim)

!-----------------------------------------------------------------------

      a13=aa**(1.0d0/3.0d0)
      rv=r0r*a13
      ri=r0i*a13
      rd=r0d*a13
      rrs=r0rs*a13
      ris=r0is*a13

      do ir=1,nrmax
       r=(ir-1)*dr
       rdv=r
       if(rdv==0.0d0) rdv=1.0d-10
       argv=(r-rv)/ar
       argi=(r-ri)/ai
       argd=(r-rd)/ad
       argrs=(r-rrs)/ars
       argis=(r-ris)/ais

       potv(ir) =-v0/(1.0d0+exp(argv))
       potw(ir) =-w0/(1.0d0+exp(argi))
       potd(ir) =-4.0d0*wd*exp(argd)/(1.0d0+exp(argd))**2
       potrs(ir)=-vso*2.0d0/ars*exp(argrs)/(1.0d0+exp(argrs))**2/rdv
       potis(ir)=-wso*2.0d0/ais*exp(argis)/(1.0d0+exp(argis))**2/rdv

      end do

      return

end subroutine wspotgen
!=======================================================================



!=======================================================================
subroutine kdpot(izp,za,aa,elab,v0,r0r,ar,w0,r0i,ai,wd,r0d,ad &
     &          ,vso,r0rs,ars,wso,r0is,ais)
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      tn=aa-za
      alp=(tn-za)/aa

      if(izp==1) go to 1000

!-----neutron depth parameters

      v1n=59.30d0-21.0d0*alp-0.024d0*aa
      v2n=0.007228d0-1.48d-6*aa
      v3n=1.994d-5-2.0d-8*aa
      v4n=7.0d-9

      w1n=12.195d0+0.0167d0*aa
      w2n=73.55d0+0.0795d0*aa

      d1n=16.0d0-16.0d0*alp
      d2n=0.018d0+0.003802d0/(1.0d0+exp((aa-156.0d0)/8.0d0))
      d3n=11.5d0

      vs1n=5.922d0+0.0030d0*aa
      vs2n=0.004d0

      ws1n=-3.1d0
      ws2n=160.0d0

      enf=-11.2814d0+0.02646d0*aa

      edf=(elab-enf)
      v0=v1n*(1.0d0-v2n*edf+v3n*edf**2-v4n*edf**3)
      w0=w1n*edf**2/(edf**2+w2n**2)
      wd=d1n*edf**2/(edf**2+d3n**2)*exp(-d2n*edf)
      vso=vs1n*exp(-vs2n*edf)
      wso=ws1n*edf**2/(edf**2+ws2n**2)

      go to 2000

 1000 continue

!-----proton depth parameters

      v1p=59.30d0+21.0d0*alp-0.024d0*aa
      v2p=0.007067d0+4.23d-6*aa
      v3p=1.729d-5+1.136d-8*aa
      v4p=7.0d-9

      w1p=14.667d0+0.009629d0*aa
      w2p=73.55d0+0.0795d0*aa

      d1p=16.0d0+16.0d0*alp
      d2p=0.018d0+0.003802d0/(1.0d0+exp((aa-156.0d0)/8.0d0))
      d3p=11.5d0

      vs1p=5.922d0+0.0030d0*aa
      vs2p=0.004d0

      ws1p=-3.1d0
      ws2p=160.0d0

      epf=-8.4075d0+0.01378d0*aa

      r0c=1.198d0+0.697d0*aa**(-2.0d0/3.0d0)+12.994d0*aa**(-5.0d0/3.0d0)
      vcb=1.73d0/r0c*za*aa**(-1.0d0/3.0d0)

      edf=(elab-epf)
      vcc=vcb*(v2p-2.0d0*v3p*edf+3.0d0*v4p*edf**2)*v1p
      v0=v1p*(1.0d0-v2p*edf+v3p*edf**2-v4p*edf**3)+vcc
      w0=w1p*edf**2/(edf**2+w2p**2)
      wd=d1p*edf**2/(edf**2+d3p**2)*exp(-d2p*edf)
      vso=vs1p*exp(-vs2p*edf)
      wso=ws1p*edf**2/(edf**2+ws2p**2)

 2000 continue

!-----geometrical parameters

      r0r=1.3039d0-0.4054d0*aa**(-1.0d0/3.0d0)
      ar=0.6778d0-1.487d-4*aa

      r0i=r0r
      ai=ar

      r0d=1.3424d0-0.01585d0*aa**(1.0d0/3.0d0)
      ad=0.5446d0-1.656d-4*aa
      if(izp==1) ad=0.5187d0+5.205d-4*aa

      r0rs=1.1854d0-0.647d0*aa**(-1.0d0/3.0d0)
      ars=0.59d0

      r0is=r0rs
      ais=ars

      return

end subroutine kdpot
!=======================================================================



!=======================================================================
!**
!**        ***** suphod *****
!**
!    this program is o.k. only when jisu=3
!                                    and dx=const  -- care --

function suphod(xin,x,dx,y,nmax,ndim)
      use kibmod
      implicit real*8(a-h,o-z)
      real*8 x(ndim),y(ndim)
      integer,parameter :: jisu=3
!** ===========================================
      if(nmax<=jisu) then
        write(kibout,400) nmax,jisu
  400   format(/1x,' nmax,jisu=',2i5,10x,'stop in suphod')
        stop 1
      end if

      nnr=int((xin-x(1))/ dx +1.01)

      if(nnr<=0 .or. nnr>nmax) then
        write(kibout,500) nnr
  500   format(/1x,' nnr=',i5,10x,'stop in suphod --- xin < x(1) ' &
     &            ,' or  xin > x(nmax) ---')
        stop 1
      end if

      jmin=nnr-1
      jmax=nnr+2
      if(jmin<=1) jmin=1
      if(jmax>=nmax) jmin=nmax-3

      x0=x(jmin)
      x1=x(jmin+1)
      x2=x(jmin+2)
      x3=x(jmin+3)

      y0=y(jmin)
      y1=y(jmin+1)
      y2=y(jmin+2)
      y3=y(jmin+3)

      xin0=xin-x0
      xin1=xin-x1
      xin2=xin-x2
      xin3=xin-x3

      x01=x0-x1
      x02=x0-x2
      x03=x0-x3

      x10=x1-x0
      x12=x1-x2
      x13=x1-x3

      x20=x2-x0
      x21=x2-x1
      x23=x2-x3

      x30=x3-x0
      x31=x3-x1
      x32=x3-x2

      xx0= xin1/x01 *xin2/x02 *xin3/x03
      xx1= xin0/x10 *xin2/x12 *xin3/x13
      xx2= xin0/x20 *xin1/x21 *xin3/x23
      xx3= xin0/x30 *xin1/x31 *xin2/x32

      suphod=xx0*y0 +xx1*y1 +xx2*y2 +xx3*y3

      return

end function suphod
!=======================================================================



!=======================================================================
!**
!**        ***** suphodx2 *****
!**

function suphodx2(xin,x,y,n1max,n2cal,n1dim,n2dim)
      use kibmod
      implicit real*8(a-h,o-z)
      real*8 x(n1dim),y(n1dim,n2dim)
!** ===========================================
      if(n1max<4) then
        write(kibout,400) n1max
  400   format(/1x,' n1max=',i5,'< 4',10x,'stop in suphodx2')
        stop 1
      end if

      if(xin<x(1)) then
        write(*,*) 'ERROR: xin<x(1)'
        write(*,'(a,f0.4,a,f0.4)') ' xin=',xin,', x(1)=',x(1)
        stop 1
      end if
      if(xin>x(n1max)) then
        write(*,*) 'ERROR: xin>x(n1max)'
        write(*,'(a,f0.4,a,f0.4)') ' xin=',xin,', x(n1max)=',x(n1max)
        stop 1
      end if

      do i=2,n1max-2
       if(xin<x(i)) then
         ibs=i-2
         if(ibs<=0) ibs=1
         go to 100
       end if
      end do
      ibs=n1max-3
  100 continue

      jmin=ibs

      x0=x(jmin)
      x1=x(jmin+1)
      x2=x(jmin+2)
      x3=x(jmin+3)

      y0=y(jmin  ,n2cal)
      y1=y(jmin+1,n2cal)
      y2=y(jmin+2,n2cal)
      y3=y(jmin+3,n2cal)

      xin0=xin-x0
      xin1=xin-x1
      xin2=xin-x2
      xin3=xin-x3

      x01=x0-x1
      x02=x0-x2
      x03=x0-x3

      x10=x1-x0
      x12=x1-x2
      x13=x1-x3

      x20=x2-x0
      x21=x2-x1
      x23=x2-x3

      x30=x3-x0
      x31=x3-x1
      x32=x3-x2

      xx0= xin1/x01 *xin2/x02 *xin3/x03
      xx1= xin0/x10 *xin2/x12 *xin3/x13
      xx2= xin0/x20 *xin1/x21 *xin3/x23
      xx3= xin0/x30 *xin1/x31 *xin2/x32

      suphodx2=xx0*y0 +xx1*y1 +xx2*y2 +xx3*y3

      return

end function suphodx2
!=======================================================================



!=======================================================================
!-----------------------------------------------------------------------
!---  subroutine fopen by Y. Iseri (Chiba-Keizai college)  ---
!-----------------------------------------------------------------------
!**
!**        ***** fopen *****
!**
subroutine fopen
      implicit real*8(a-h,o-z)
      character fname*50,sta*8,comm*60,off*1
! ==========================================
! --- unit-5 ---
!     write(*,100)
! 100 format(//1x,'>>> input file-name of unit 5 ==> ',$)
!     read(*,'(a)') fname

!     open(5,file=fname,status='old')

! --- other units ---
      read(5,'(a)') comm
      write(*,200) comm
  200 format(/1x,'---- (comment in file) ----' &
     &       /1x,' << ',a,' >>'                &
     &      //1x,'---- (  open  files  ) ----' &
     &      //1x,'  unit  staus      file-name')

    1 read(5,300) off,iunit,sta,fname
  300 format(a1,i3,1x,a8,2x,a50)
       if(off/=' ') go to 1
       if(iunit>=100 .or. iunit<0) go to 900
       write(*,310) iunit,sta,fname
  310  format(1x,3x,i2,3x,a,3x,a)
       open(iunit,file=fname,status=sta)
      go to 1

! --- one more line for discrimination ---
  900 read(5,'(a)') comm
      return

end subroutine fopen
!=======================================================================



!=======================================================================
!-----------------------------------------------------------------------
!---  subroutine element by Y. Iseri (Chiba-Keizai college)  ---
!-----------------------------------------------------------------------
!**
!**        ***** elemnt *****
!**
subroutine elemnt(amass,z,mass,name)
      implicit real*8(a-h,o-z)
      integer,parameter :: mmax=109
      character*2 name
      character*2 nat(mmax),nap(3),naa(2),neu,non
      data  neu,nat / 'n ',                                         &
     &      'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
     &      'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti', &
     &      'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As', &
     &      'Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru', &
     &      'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs', &
     &      'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
     &      'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
     &      'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
     &      'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es', &
     &      'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt'/
      data nap /'p ','d ','t '/,  naa/'-p','-n'/, non/'??'/
!** ============================================================
      mass=nint(amass)
      iz=nint(z)

      if(iz==-1) then
         name= naa(1)
      else if(iz== 0) then
         name= neu
      else if(iz== 1 .and. mass<=3) then
         name= nap(mass)
      else if(iz>=1 .and. iz<=mmax) then
         name= nat(iz)
      else
         name= non
      end if

      return

end subroutine elemnt
!=======================================================================



!=======================================================================
subroutine tdxqm(rmax,dr,fkayi,thtkg,phikg                    &
     &          ,elab,amasspi,e0,e1,e2,ea,eb,t0a,t1a,t2a,ylm0 &
     &          ,zp,z1,z2,za,zb,pwsum3                        &
     &          ,rcn,ircn,pmi,beta,sfac,ucou,dthnn,dsigelm    &
     &          ,wf,wfmtp,ffrg,nrmax,i1st,fj,ls,lmax,icpot    &
     &          ,ngth,ngph,ielm,kinelm,sig                    &
     &          ,ipot,facv,facw,facvs,facws,eelm,iex,io       &
     &          ,izsp,ap,asp,aa,ab,th,costh                   &
     &          ,fnl,fnltbl,ims,iedg,ielmedg                  &
     &          ,iepfix,epfix,nepfmax,epfmin                  &
     &          ,detbl,uopttbl,ulstbl,uopt                    &
     &          ,potv,potw,potd,potrs,potis                   &
     &          ,fm0,fm1,fm2,fmsp                             &
     &          ,nthnnmax,thnnmin,nennmax,sigelm              &
     &          ,fk0,th0,ph0,fk1,th1,ph1,fk2,th2,ph2          &
     &          ,thcalg,thweig,nthgmax                        &
     &          ,phcalg,phweig,nphgmax                        &
     &          ,rcalg,rweig,nrgmax                           &
     &          ,pltbl,dpltbl,dzpl,nzplmax                    &
     &          ,hv0,c0r,fmuelm,const,tdxr                    &
     &          ,itmdcal,ilscal,wfls,uls,ay)
      use consts
      use dims
      use array, only: gma
      use kibmod
      use nntbl
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 fkayi(3),thtkg(3),phikg(3),amasspi(3)
      real*8 ylm0(nthdim,ldim,mdim),pmi(3),beta(3)
      real*8 ucou(3,nrdim),ffrg(ngdim)
      real*8 dsigelm(nenndim,nthnndim),thelm(nthnndim)
      real*8 facv(3),facw(3),facvs(3),facws(3)
      integer ims(3),iedg(3)
      real*8 eelm(nenndim)

      complex*16 wf(3,ngdim,ngdim,ngdim)
      complex*16 wfmtp(ngdim,ngdim,ngdim)

      integer lmax(3),icpot(3),ipot(3)
      integer iepfix(3),nepfmax(3)
      real*8 epfix(3),epfmin(3)
      real*8 th(nthdim),costh(nthdim)
      real*8 potv(nrdim),potw(nrdim),potd(nrdim),potrs(nrdim),potis(nrdim)
      complex*16 uopt(3,nrdim),uopttbl(3,nepfdim,nrdim)
      complex*16 ulstbl(3,nepfdim,nrdim)
      complex*16 uoptls(3,nrdim)
      complex*16 wfls(2,2,3,ngdim,ngdim,ngdim)
      complex*16 uls(3,nrdim)
      complex*16 allsumm(-mdim:mdim)
      complex*16 tmtxls(2,2,2,2,2,2,-mdim:mdim)
      complex*16 ylmphi2,ylmphi2tbl(-ldim:ldim,ngdim)
      real*8 sigy(2)
      complex*16 elmnntbl(16)
      complex*16 fnl(3,nrdim),fnltbl(3,nepfdim,nrdim)

      real*8 thcalg(ngdim),thweig(ngdim),rcalg(ngdim)
      real*8 phcalg(ngdim),phweig(ngdim),rweig(ngdim)
      real*8 pltbl(ldim,nzpldim),dpltbl(ldim,nzpldim)

      character*1 :: excd(6)=' '

!-----------------------------------------------------------------------

      tdxr=0
      ay=0.d0
      if (fk1<0.02d0) return
      if (fk2<0.02d0) return

      if(i1st==1) then


!----------------------------------------------------
!---  preparation for Gauss Legendre integration  ---
!----------------------------------------------------

        nthgmax=ngth
        call glweight(nthgmax,0.d0,pi,thcalg,thweig)

        nphgmax=ngph
        call glweight(nphgmax,0.d0,2.d0*pi,phcalg,phweig)


!----------------------------------------------
!---  making a table of plfunc and dplfunc  ---
!----------------------------------------------

        lmaxmax=maxval(lmax(:))
        dzpl=1.0d-3
        nzplmax=nint(2.0d0/dzpl)+1
        if(nzplmax>nzpldim) then
          write(*,*) 'ERROR: nzplmax > nzpldim'
          write(*,'(a,i0,a,i0)') ' nzplmax=',nzplmax, ', nzpldim=',nzpldim
          stop 1
        end if

        do l=0,lmaxmax
          l1=l+1
          do iz=1,nzplmax
           z=(iz-1)*dzpl-1.0d0
           pltbl(l1,iz)=plfunc(l,z)
           dpltbl(l1,iz)=dplfunc(l,z)
          end do
         end do

        rcn=rmax-3.0d0*dr
        ircn=nint(rcn/dr)+1
        pmi(1)=1.0d0
        pmi(2)=-1.0d0
        pmi(3)=-1.0d0

        do j=1,3
         if(ipot(j)==0) then
           icpot(j)=0
           facv(j)=0.0d0
           facw(j)=0.0d0
         end if

        end do

        amasspi(1)=ap
        amasspi(2)=ap
        amasspi(3)=asp
        if(iex==1) then
          amasspi(2)=asp
          amasspi(3)=ap
        end if

!--for Gauss-Legendre nodes
        lspmax=ls
        mspmax=ls
        nthmin=1
        nthmax=nthgmax
        do ith=nthmin,nthmax
         th(ith)=thcalg(ith)
         rad=th(ith)
         costh(ith)=cos(rad)
        end do
!--for Gauss-Legendre nodes end

        call ylmph0(ldim,mdim,nthdim,lspmax,mspmax,nthmin,nthmax,costh,ylm0)

        fkayi(1)=fk0
        thtkg(1)=th0*piad
        phikg(1)=ph0*piad


!----------------------------------
!---  potential for particle 0  ---
!----------------------------------

        izp=nint(zp)
        uoptls(1,:)=uls(1,:)
        if(ipot(1)>9) then
          call potint(nrmax,nepfmax(1),iepfix(1),iedg(1)   &
     &               ,epfmin(1),epfix(1),t0a,detbl,beta(1) &
     &               ,uopttbl(1,1:nepfmax(1),1:nrmax)      &
     &               ,ulstbl(1,1:nepfmax(1),1:nrmax)       &
     &               ,fnltbl(1,1:nepfmax(1),1:nrmax)       &
     &               ,uopt(1,1:nrmax),uoptls(1,1:nrmax),fnl(1,1:nrmax))
         else if(ipot(1)==1) then
          call kdpot(izp,za,aa,elab,v0,r0r,ar,w0,r0i,ai &
     &              ,wd,r0d,ad,vso,r0rs,ars,wso,r0is,ais)
          call wspotgen(aa,v0,r0r,ar,w0,r0i,ai,wd,r0d,ad &
     &                 ,vso,r0rs,ars,wso,r0is,ais        &
     &                 ,potv,potw,potd,potrs,potis,nrmax,dr)
          do ir=1,nrmax
           uopt(1,ir)=facv(1)*potv(ir)+wi*facw(1)*(potw(ir)+potd(ir))
           uoptls(1,ir)=facvs(1)*potrs(ir)+wi*facws(1)*potis(ir)
          end do
        end if

        if(kinelm==1) then
          if(ielm/=4) then
            write(kibout,601)
  601       format(/3x,'-- kinematics of the elementary process --' &
     &            //8x,'config#',3x,'fkapf', 5x,'fkapi', 5x,'thnn'  &
     &             ,6x,'fkapfe',4x,'fkapie',4x,'thnne'              &
     &             ,4x,'ennin',5x,'beteg',5x,'fmol2')
           else
            write(kibout,602)
  602       format(/3x,'-- kinematics of the elementary process --' &
     &            //8x,'config#',3x,'fkapf', 5x,'fkapi', 5x,'thnn'  &
     &             ,5x,'fkapfe',4x,'fkapie',5x,'thnne'              &
     &             ,4x,'ennin',5x,'beteg',5x,'fmol2'                &
     &             ,4x,'thkapfe',3x,'phkapfe',3x,'thkapie',3x,'phkapie')
          end if
        end if

        call sigelm0(ielm,elab,izp,izsp,ap,asp,nthnnmax,thnnmin,dthnn,nennmax &
     &              ,sigelm,eelm,dsigelm)

        hv0=hc**2*fk0*(e0+ea)/(e0*ea+hc**2*fk0**2)
        c0r=1.0d0/(hv0*(2.0d0*pi)**5*(2.0d0*ls+1.0d0))
        fmuelm=fm1*fm2/(fm1+fm2)

        if(ielm==0) then
          sigfac=sigelm
         else if(ielm==1 .or. ielm==3) then
          sigfac=1.d0
         else if(ielm==4) then
          sigfac=1.d0
          c0r=1.0d0/(hv0*(2.0d0*pi)**5*2.0d0*(2.0d0*fj+1.0d0))
         else
          sigfac=1.0d0
        end if

        const=c0r*sigfac*sfac*10.0d0  ! fm^2 --> mb

      end if ! if(i1st==1) end


!----------------------------------------------------------
!---  kinematics and table for the elementary process  ---
!----------------------------------------------------------

      if((ielm==3 .or. ielm==4) .or. kinelm==1)             &
     &  call kinemaelm(1,fk1,e1,th1,ph1,fk2,e2,th2,ph2,aa   &
     &                ,fk0,th0,ph0,fm0,fmsp,fm1,fm2,iex     &
     &                ,fkapf,fkapi,thnn,fkapfe,fkapie,thnne &
     &                ,thkapie,thkapfe,phkapie,phkapfe      &
     &                ,ennin,beteg,fmol2,fkaponsh,fmurel)

      if(kinelm==1 .and. ielm/=4)                                 &
     &  write(kibout,603) io,fkapf,fkapi,thnn,fkapfe,fkapie,thnne &
     &              ,ennin,beteg,fmol2
  603 format(4x,i10,7f10.4,2f10.5)

      if(iex==1) thnn =180.0d0-thnn
      if(iex==1) thnne=180.0d0-thnne

      fkayi(2)=fk1
      fkayi(3)=fk2
      thtkg(2)=th1*piad
      thtkg(3)=th2*piad
      phikg(2)=ph1*piad
      phikg(3)=ph2*piad

      nq0cal  =0
      nqcal   =0
      nthq0cal=0
      nthqcal =0
      nphq0cal=0
      nphqcal =0

      if(ielm==4) then
        if(ifnntbl==1) then
          read(kibelm,*) q0mn,q0mx,dq0
          read(kibelm,*) th0mn,th0mx,thmn,thmx,dth
          ph0mn=0.d0
          ph0mx=360.d0
          phmn=0.d0
          phmx=360.d0
          dph=180.d0

          if(q0mn>=q0mx) then
            write(*,*) 'ERROR: q0mn must be smaller than q0mx'
            write(*,'(a,f0.4,a,f0.4)') ' q0mn',q0mn,', q0mx',q0mx
            stop 1
          end if
          if(dq0<0.0d0) then
            write(*,*) 'ERROR: dq0 must be positive'
            write(*,'(a,f0.4)') ' dq0=',dq0
            stop 1
          end if
          if(th0mn>=th0mx) then
            write(*,*) 'ERROR: th0mn must be smaller than th0mx'
            write(*,'(a,f0.4,a,f0.4)') ' th0mn',th0mn,', th0mx',th0mx
            stop 1
          end if
          if(thmn>=thmx) then
            write(*,*) 'ERROR: thmn must be smaller than thmx'
            write(*,'(a,f0.4,a,f0.4)') ' thmn',thmn,', thmx',thmx
            stop 1
          end if
          if(dth<0.0d0) then
            write(*,*) 'ERROR: dth must be positive'
            write(*,'(a,f0.4)') ' dth=',dth
            stop 1
          end if
          if(ph0mn>=ph0mx) then
            write(*,*) 'ERROR: ph0mn must be smaller than ph0mx'
            write(*,'(a,f0.4,a,f0.4)') ' ph0mn',ph0mn,', ph0mx',ph0mx
            stop 1
          end if
          if(phmn>=phmx) then
            write(*,*) 'ERROR: phmn must be smaller than phmx'
            write(*,'(a,f0.4,a,f0.4)') ' phmn',phmn,', phmx',phmx
            stop 1
          end if
          if(dph<0.0d0) then
            write(*,*) 'ERROR: dph must be positive'
            write(*,'(a,f0.4)') ' dph=',dph
            stop 1
          end if

          qmn=q0mn
          qmx=q0mx
          dq=dq0

          call makenntbl()
          ! All arguments to be kept when (ifnntbl/=1) are in module 'nntbl'.

          ! read in NN input and make on-shell NN t-matrix table
          call onshnnread(nthqdim,nphqdim                         &
     &                   ,nq0dim,nthq0dim,nphq0dim,gma            &
     &                   ,nq0mn,nq0mx,nth0mn,nth0mx,nph0mn,nph0mx &
     &                   ,nthmn,nthmx,nphmn,nphmx)
          ifnntbl=0
        end if

!------- on shell approx
        q0=fkaponsh ! fkaponsh is calculated in kinemaelm
        q=q0

        thq0=thkapie*piad
        thq=thkapfe*piad
        phq0=phkapie*piad
        phq=phkapfe*piad

        call nnkincal(q0,q,thq0,thq,phq0,phq,ielmedg,ionsh ,excd &
     &               ,nq0cal,nqcal,nphq0cal,nphqcal              &
     &               ,q0cal,qcal,thq0cal,thqcal,phq0cal,phqcal)

        if(q<qmn .or. q>qmx) then
          excd(1)='%'
          if(ionsh==0 .or. ionsh==1) excd(1)='*'
        end if
        if(q0<q0mn .or. q0>q0mx) then
          excd(2)='%'
          if(ionsh==0 .or. ionsh==2) excd(2)='*'
        end if

        if(kinelm==1)                                              &
     &    write(kibout,604) io,fkapf,fkapi,thnn,fkapfe,excd(1)     &
     &                     ,fkapie,excd(2),thnne,ennin,beteg,fmol2 &
     &                     ,thq/piad,excd(3),phq/piad,excd(4)      &
     &                     ,thq0/piad,excd(5),phq0/piad,excd(6)
  604   format(4x,i10,3f10.4,2(f9.4,a1),2f10.4,2f10.5,4(f9.4,a1))

      end if ! if(ielm==4) end


!---------------------------------------
!---  potential for particle 1 and 2 ---
!---------------------------------------

      t1=t1a
      t2=t2a
      iz1=nint(z1)
      iz2=nint(z2)
      do idw=2,3
       if(idw==2) then
         tin=t1
         izin=iz1
        else if(idw==3) then
         tin=t2
         izin=iz2
       end if
       uoptls(idw,:)=uls(idw,:)

       if(ipot(idw)>9) then
         call potint(nrmax,nepfmax(idw),iepfix(idw),iedg(idw)     &
     &                ,epfmin(idw),epfix(idw),tin,detbl,beta(idw) &
     &                ,uopttbl(idw,1:nepfmax(idw),1:nrmax)        &
     &                ,ulstbl(idw,1:nepfmax(idw),1:nrmax)         &
     &                ,fnltbl(idw,1:nepfmax(idw),1:nrmax)         &
     &                ,uopt(idw,1:nrmax),uoptls(idw,1:nrmax)      &
     &                ,fnl(idw,1:nrmax))
        else if(ipot(idw)==1) then
         call kdpot(izin,zb,ab,tin,v0,r0r,ar,w0,r0i,ai,wd,r0d,ad &
     &             ,vso,r0rs,ars,wso,r0is,ais)
         call wspotgen(ab,v0,r0r,ar,w0,r0i,ai,wd,r0d,ad &
     &                ,vso,r0rs,ars,wso,r0is,ais        &
     &                ,potv,potw,potd,potrs,potis,nrmax,dr)

         do ir=1,nrmax
          uopt(idw,ir)=facv(idw)*potv(ir)+wi*facw(idw)*(potw(ir)+potd(ir))
          uoptls(idw,ir)=facvs(idw)*potrs(ir)+wi*facws(idw)*potis(ir)
         end do

       end if
      end do

      do idw=1,3
       if(ipot(idw)==0) then
         uopt(idw,:)=(0.0d0,0.0d0)
         ucou(idw,:)=0.0d0
         fnl(idw,:)=1.0d0
         uoptls(idw,:)=(0.0d0,0.0d0)
       end if
      end do

      call pwtmtx(nrgmax,rcalg,rweig,ls,fkayi,thtkg,phikg,amasspi,aa,pwsum3,ffrg,iex)
      if(ielm /= 4 .and. ipot(1)+ipot(2)+ipot(3) == 0) then
        allmsum=pwsum3
      else
        wfmtp(1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
        if(i1st==1) then  ! first time calc. of t-matrix
          idwmin=1
          wf(:,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
          wfls(:,:,:,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
         else
          idwmin=2
          wf(2,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
          wf(3,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
          wfls(:,:,2,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
          wfls(:,:,3,1:nrgmax,1:nthgmax,1:nphgmax)=(0.0d0,0.0d0)
        end if

        do idw=idwmin,3

         amassp=amasspi(idw)

         if(idw==1) then
           amasst=aa
           zpt=zp
           ztg=za
           rede=e0*ea/(e0+ea)
          else if(idw==2) then
           amasst=ab
           zpt=z1
           ztg=zb
           rede=e1*eb/(e1+eb)
          else if(idw==3) then
           amasst=ab
           zpt=z2
           ztg=zb
           rede=e2*eb/(e2+eb)
          else
           write(*,*) 'error in idw selection'
           stop 1
         end if

         if(ims(idw)==0) then
           fmu=rede/ac
          else
           fmu=(amassp*amasst)/(amassp+amasst)
         end if

         call dwcalc(ircn,nrmax,ilscal,nthgmax,nphgmax,nrgmax &
     &              ,icpot(idw),lmax(idw),nzplmax,idw         &
     &              ,dr,rcn,zpt,ztg,fmu,fkayi(idw),thtkg(idw) &
     &              ,phikg(idw),pmi(idw),beta(idw),dzpl       &
     &              ,thcalg,phcalg,rcalg,ucou(idw,:),pltbl    &
     &              ,uopt(idw,:),uoptls(idw,:),fnl(idw,:)     &
     &              ,wfls(:,:,idw,:,:,:)  )

        end do

        icalls1=2
        icalls2=2
        icalls3=2
        if(ilscal==0) then
          icalls1=1
          icalls2=1
          icalls3=1
        end if
        allmsumls=0.0d0
        tmtxls(:,:,:,:,:,:,:)=(0.0d0,0.0d0)
        ncal=0
        mvalmin=-ls
        if((abs(phikg(2))<epslim .and. abs(phikg(3)-pi)<epslim) .or. &
     &     (abs(phikg(2)-pi)<epslim .and. abs(phikg(3))<epslim)) mvalmin=0
        if(ilscal==1 .or. ielm==4) mvalmin=-ls

        do mval=-ls,ls
         do iphg=1,nphgmax
          ylmphi2tbl(mval,iphg)=ylmphi2(mval,phcalg(iphg))
         end do
        end do

        do iup1i=1,icalls1
         do iup2i=1,icalls2
          do iup3i=1,icalls3
           do iup1f=1,icalls1
            if(iup1f==2 .and. iup1i==2) cycle ! non-flip compontents are the same
            do iup2f=1,icalls2
             if(iup2f==2 .and. iup2i==2) cycle ! non-flip compontents are the same
             do iup3f=1,icalls3
              if(iup3f==2 .and. iup3i==2) cycle ! non-flip compontents are the same

              do ir=1,nrgmax
               do ithg=1,nthgmax
                do iphg=1,nphgmax

                 wfmtp(ir,ithg,iphg)=wfls(iup1i,iup1f,1,ir,ithg,iphg) &
     &                              *wfls(iup2i,iup2f,2,ir,ithg,iphg) &
     &                              *wfls(iup3i,iup3f,3,ir,ithg,iphg)
                end do
               end do
              end do
!---t-matrix calc. for each m0, m'0, m1, m'1, m2 and m2'.
              allsumm=(0.0d0,0.0d0)
              call tmtx(dr,fkayi,thtkg,amasspi,aa,ylm0,wfmtp               &
     &                 ,ffrg,nrmax,ls,thcalg,thweig,nthgmax,phweig,nphgmax &
     &                 ,rcalg,rweig,nrgmax,iex,io,itmdcal,allsumm,mvalmin  &
     &                 ,ylmphi2tbl)

              do mval=mvalmin,ls
               tmtxls(iup1i,iup1f,iup2i,iup2f,iup3i,iup3f,mval)=allsumm(mval)
              end do
              ncal=ncal+1

             end do ! up3f end
            end do ! iup2f end
           end do ! iup1f end
          end do ! iup3i end
         end do ! iup2i end
        end do ! iup1i end

!-------- non-flip compontents are the same ---
        tmtxls(2,2,:,:,:,:,:)=tmtxls(1,1,:,:,:,:,:)
        tmtxls(:,:,2,2,:,:,:)=tmtxls(:,:,1,1,:,:,:)
        tmtxls(:,:,:,:,2,2,:)=tmtxls(:,:,:,:,1,1,:)


        sallsum=0.d0

!-------- NN cross section approximation ---
        if (ielm/=4) then
          do mval=mvalmin,ls
           fmfac=1.0d0
           if(mvalmin==0 .and. mval/=0) fmfac=2.0d0
           sallsum=sallsum+fmfac*abs(tmtxls(1,1,1,1,1,1,mval))**2.0d0
          end do

          ay=0.d0 ! Ay=0.0 is returned when the NN c.s. approx. is adopted

!-------- NN t-matrix ---
         else if(ielm==4) then
          call getelmnntbl(nphqcal,nq0cal,nphq0cal,thqcal,thq0cal,elmnntbl)

          call sigycalc(ls,thtkg,phikg,fj,elmnntbl,tmtxls,sigy)

          sallsum=sum(sigy(:))
          ay=(sigy(1)-sigy(2))/sallsum
        end if

        faclscal=dble(icalls1*icalls2*icalls3)

        if(ielm==4) faclscal=1.0d0
        allmsum=sallsum/faclscal

      end if

      if(ielm==0) then
        sig=1.0d0
       else if(ielm==3) then
        nthnne=nint(thnne/dthnn)+1
        enninc=ennin
        if(ennin<eelm(1)) then
          if(ielmedg==0) then
            write(*,*) 'ERROR: ennin<eelm(1)'
            write(*,'(a,f0.4,a,f0.4)') ' ennin=',ennin,' eelm(1)=',eelm(1)
            stop 1
          end if
          enninc=eelm(1)
         else if(ennin>eelm(nennmax)) then
          if(ielmedg==0) then
            write(*,*) 'ERROR: ennin>eelm(nennmax)'
            write(*,'(a,f0.4,a,f0.4)') ' ennin=',ennin,' eelm(nennmax)=',eelm(nennmax)
            stop 1
          end if
          enninc=eelm(nennmax)
        end if
        do ithelm=1,nthnnmax
          thelm(ithelm)=(ithelm-1)*dthnn+thnnmin
        end do
        sig=polintm2d(eelm(1:nennmax),thelm(1:nthnnmax),dsigelm(1:nennmax,1:nthnnmax) &
     &               ,nennmax,nthnnmax,enninc,thnne,4,4)*fmol2
       else if(ielm==4) then
        sig=fmol2
      end if

      tdxr=allmsum*const*sig

      i1st=0

      return

end subroutine tdxqm
!=======================================================================



!=======================================================================
subroutine sigycalc(ls,thtkg,phikg,fj,elmnntbl,tmtxls,sigy)
      use consts,only:pi,wi
      use dims,only:mdim
      implicit none
!=======================================================================
      integer,intent(in):: ls
      real*8,intent(in):: thtkg(3),phikg(3),fj
      complex*16,intent(in):: elmnntbl(16),tmtxls(2,2,2,2,2,2,-mdim:mdim)

      real*8,intent(out):: sigy(2)

      integer iup1y,jz,jzmin,jzmax,iup2ik,iup3ik,iup1iz,iupsp,mval &
     &       ,iup1fz,iup2fk,iup3fk,ispin,iup2fz,iup3fz
      real*8 up1y,fjz,up2ik,up3ik,cgcal,dmval,up1fz,up2fk,up3fk &
     &      ,cos2kz,cos3kz,dlval,up2fz,up3fz,up1iz,upsp,clebsh,dfunc
      complex*16 sfsum

!-----------------------------------------------------------------------

      sigy(:)=0.0d0
      cos2kz=cos(thtkg(2))
      cos3kz=cos(thtkg(3))

      dlval=dble(ls)
      jzmin=nint(-2.0d0*fj)
      jzmax=nint( 2.0d0*fj)

      do iup1y=1,2
       if(iup1y==1) up1y=0.5d0
       if(iup1y==2) up1y=-0.5d0

       do jz=jzmin,jzmax,2
        fjz=jz/2.0d0

        do iup2ik=1,2
         if(iup2ik==1) up2ik=0.5d0
         if(iup2ik==2) up2ik=-0.5d0

         do iup3ik=1,2
         if(iup3ik==1) up3ik=0.5d0
         if(iup3ik==2) up3ik=-0.5d0

         sfsum=(0.0d0,0.0d0)

         do iup1iz=1,2
          if(iup1iz==1) up1iz=0.5d0
          if(iup1iz==2) up1iz=-0.5d0

          do iupsp=1,2
           if(iupsp==1) upsp=0.5d0
           if(iupsp==2) upsp=-0.5d0
           dmval=fjz-upsp
           mval=nint(dmval)
           if(mval>ls) cycle

           cgcal=clebsh(dlval,dmval,0.5d0,upsp,fj,fjz)

           do iup1fz=1,2
            if(iup1fz==1) up1fz=0.5d0
            if(iup1fz==2) up1fz=-0.5d0

            do iup2fk=1,2
             if(iup2fk==1) up2fk=0.5d0
             if(iup2fk==2) up2fk=-0.5d0

             do iup3fk=1,2
              if(iup3fk==1) up3fk=0.5d0
              if(iup3fk==2) up3fk=-0.5d0

              do iup2fz=1,2
               if(iup2fz==1) up2fz=0.5d0
               if(iup2fz==2) up2fz=-0.5d0

               do iup3fz=1,2
                if(iup3fz==1) up3fz=0.5d0
                if(iup3fz==2) up3fz=-0.5d0

                ispin=(iup1fz-1)*8+(iupsp-1)*4+(iup2fz-1)*2+(iup3fz-1)*1+1

                sfsum=sfsum+elmnntbl(ispin)*cgcal                            &
     &               *dfunc(0.5d0,up2fz,up2fk,cos2kz)*exp(wi*up2fz*phikg(2)) &
     &               *dfunc(0.5d0,up3fz,up3fk,cos3kz)*exp(wi*up3fz*phikg(3)) &
     &               *tmtxls(iup1iz,iup1fz,iup2ik,iup2fk,iup3ik,iup3fk,mval) &
     &               *dfunc(0.5d0,up1iz,up1y,0.0d0)*exp(-wi*up1iz*pi/2.d0)

                end do ! iup3fz
               end do ! iup2fz
              end do ! iup3fk
             end do ! iup2fk
            end do ! iup1fk
           end do ! iupsp
          end do ! iup1iz

          sigy(iup1y)=sigy(iup1y)+abs(sfsum)**2.0d0

         end do ! iup3ik
        end do ! iup2ik
       end do ! jz
      end do ! iup1y

end subroutine sigycalc
!=======================================================================



!=======================================================================
subroutine pwtmtx(nrgmax,rcalg,rweig,ls,fkayi,thtkg,phikg,amasspi,aa,pwsum3,ffrg,iex)
      use consts
      use dims
      implicit real*8(a-h,o-z)
!=======================================================================
      real*8 fkayi(3),thtkg(3),phikg(3)
      real*8 amasspi(3),ffrg(ngdim)
      real*8 rcalg(ngdim),rweig(ngdim)
      real*8 bj(ls+1),bn(ls+1)

!-----------------------------------------------------------------------

      thtk0=thtkg(1)
      thtk1=thtkg(2)
      thtk2=thtkg(3)
      phik0=phikg(1)
      phik1=phikg(2)
      phik2=phikg(3)
      fkay0=fkayi(1)
      fkay1=fkayi(2)
      fkay2=fkayi(3)

      qx=fkay0*sin(thtk0)*cos(phik0)-fkay1*sin(thtk1)*cos(phik1) &
     &  -fkay2*sin(thtk2)*cos(phik2)
      qy=fkay0*sin(thtk0)*sin(phik0)-fkay1*sin(thtk1)*sin(phik1) &
     &  -fkay2*sin(thtk2)*sin(phik2)
      qz=fkay0*cos(thtk0)-fkay1*cos(thtk1)-fkay2*cos(thtk2)

      if(iex==0) then
        qz=qz-fkay0*cos(thtk0)*amasspi(3)/aa
       else if(iex==1) then
        qz=qz-fkay0*cos(thtk0)*amasspi(2)/aa
      end if

      qnrm=sqrt(qx**2.0d0+qy**2.0d0+qz**2.0d0)

      pwsum=0.0d0

      do ir=1,nrgmax
       r=rcalg(ir)
       z=qnrm*r

       call coulfg(ls,0.0d0,z,1.0d0,bj(ls+1),bn(ls+1),dummy1,dummy2)
       fjl=bj(ls+1)/z

       sigr=1.0d0
       pwint=fjl*sigr*ffrg(ir)*r**2.0d0
       pwsum=pwsum+pwint*rweig(ir)

      end do

      pwsum2=abs(pwsum)**2.0d0

      pwsum3=4.0d0*pi*(2*ls+1)*pwsum2

      return

end subroutine pwtmtx
!=======================================================================



!=======================================================================
subroutine dwcalc(ircn,nrmax,ilscal,nthgmax,nphgmax,nrgmax         &
     &           ,icpot,lmaxw,nzplmax,idw                          &
     &           ,dr,rcn,zpt,ztg,fmu,fkay,thtk,phik,pm,beta,dzpl   &
     &           ,thcalg,phcalg,rcalg,ucou,pltbl,uopt,uoptls,fnl   &
     &           ,wfls)
      use consts, only: hc,ec,ac
      use dims, only: ldim,nrdim,ngdim,nzpldim
      use array,only:lauto
      implicit real*8(a-h,o-z)
!=======================================================================
      integer,intent(in) :: ircn,nrmax,ilscal,nthgmax,nphgmax,nrgmax            &
     &                     ,icpot,lmaxw,nzplmax,idw
      real*8,intent(in) :: dr,rcn,zpt,ztg,fmu,fkay,thtk,phik,pm,beta,dzpl       &
     &                    ,thcalg(ngdim),phcalg(ngdim),rcalg(ngdim),ucou(nrdim) &
     &                    ,pltbl(ldim,nzpldim)
      complex*16,intent(in) :: uopt(nrdim),uoptls(nrdim),fnl(nrdim)
      complex*16,intent(out) :: wfls(2,2,ngdim,ngdim,ngdim)

      real*8 sigma(ldim)
      complex*16 uwfsvls(ldim,nrdim,2)
      complex*16 smtxls(ldim,2)
      complex*16 childkruug(ldim,ngdim),childkrudg(ldim,ngdim)
      complex*16 pfacsvg(ngdim),exsig(ldim),pfaccal
      complex*16 hp(ldim),hm(ldim),hpd(ldim),hmd(ldim)
      complex*16 chilud(ldim,nrdim,2),chiludg(ldim,ngdim,2)

!-----------------------------------------------------------------------

      lmincal=0
      lmaxcal=lmaxw
      if(lauto(idw) == 1) then
        if(lmaxcal > nint(fkay*(nrmax-1)*dr)) then
          lmaxcal=nint(fkay*(nrmax-1)*dr)
          if(lmaxcal+1>ldim) then
            write(*,*) 'ERROR: lmax > ldim'
            write(*,'(a,i0,a,i0,a,i0)') ' lmaxcal+1=',lmaxcal+1,', idw=',idw &
     &                                 ,', ldim=',ldim
            stop 1
          end if
        end if
      end if

      ecm=(hc*fkay)**2.0d0/(2.0d0*fmu*ac)
      fm=2.0d0*fmu*ac/hc**2

!--- Coulomb phase shift ---
      eta=zpt*ztg*fmu*ac/(ec*hc*fkay)
      rho=fkay*rcn
      drho=fkay*dr
      if(icpot==1) then
        call sigmal(eta,lmaxcal+1,sigma,exsig)
       else
        eta=0.d0
        sigma(:)=0.0d0
      end if


!-----------------------
!---  wave function  ---
!-----------------------

!--- Radial w.f. calculation (unnormalized) ---
      call wfqm(lmaxw,nrmax,ilscal,icpot,ecm,fm,dr,ucou,uopt,uoptls,uwfsvls)

!--- Wave function at the connecting point ---
      do l=1,lmaxcal+1
       call coulfg(l-1,eta,rho,drho,f1,g1,fp1,gp1)
       hp(l)=dcmplx(g1,f1)
       hm(l)=dcmplx(g1,-f1)
       hpd(l)=dcmplx(gp1,fp1)*fkay
       hmd(l)=dcmplx(gp1,-fp1)*fkay
      end do

!--- S-matrix calculation ---
      do i=1,2
       call deltacalc(lmaxcal,ircn,dr,hp,hm,hpd,hmd,uwfsvls(:,:,i),smtxls(:,i))
      end do

!--- Radial w.f. normalization ---
      call wfnrm(lmaxcal,nrmax,ircn,fkay,dr,hp,hm,smtxls,uwfsvls,chilud)

!--- Interpolation of radial w.f. at the Legendre nodes ---
      call leg_inpl(nrmax,nrgmax,lmaxcal,dr,rcalg,chilud,chiludg)

!--- making non-flip and flipped radial wave functions ---
      call makechil(lmaxcal,nrgmax,fkay,sigma,rcalg,chiludg,childkruug,childkrudg)

!--- Making 3-dimensional distorted wave ---
      call wfout(lmincal,lmaxcal,nthgmax,nphgmax,nrgmax,nzplmax,ilscal &
     &          ,fkay,thtk,phik,pm,thcalg,phcalg,rcalg,pltbl,dzpl      &
     &          ,wfls,childkruug,childkrudg)

!--- Non-local correction for the distorted wave ---
      if(abs(beta)>0.d0) then
        call make_pfacsvg(nrmax,nrgmax,beta,dr,fmu,uopt,fnl,rcalg,pfacsvg)
        do ir=1,nrgmax
         pfaccal = pfacsvg(ir)
         if(pm < 0.d0) then 
           pfaccal = conjg(pfaccal)
         end if
         wfls(:,:,ir,1:nthgmax,1:nphgmax)=wfls(:,:,ir,1:nthgmax,1:nphgmax)*pfaccal
        end do
      end if

      return

end subroutine dwcalc
!=======================================================================



!=======================================================================
subroutine wfqm(lmaxw,nrmax,ilscal,icpot,ecm,fm,dr,ucou,uopt,uoptls,uwfsvls)
      use consts,only:hc,ac,epslim
      use dims,only:nrdim,ldim
      use kibmod,only:kibout
      implicit real*8 (a-h,o-z)
!=======================================================================
      integer,intent(in) :: lmaxw,nrmax,ilscal,icpot
      real*8,intent(in) :: ecm,fm,dr,ucou(nrdim)
      complex*16,intent(in) :: uopt(nrdim),uoptls(nrdim)

      complex*16,intent(out) :: uwfsvls(ldim,nrdim,2)

      complex*16 uoptcal(nrdim)
      complex*16 optin(nrdim),fnumer(nrdim)

!-----------------------------------------------------------------------

      if(icpot==1) then
        uoptcal(1:nrmax)=uopt(1:nrmax)+dcmplx(ucou(1:nrmax),0.d0)
       else
        uoptcal(1:nrmax)=uopt(1:nrmax)
      end if

      uwfsvls(:,:,:)=(0.0d0,0.0d0)
      flsfac=0.d0

      do ijud=1,ilscal+1
       do ll=0,lmaxw
        fll1=dble(ll*(ll+1))

        if(ilscal==1) then
          if(ijud==1) then
            lsfac=ll
           else if(ijud==2) then
            lsfac=-(ll+1)
           else
            write(kibout,*) 'error in ijud'
            stop 1
          end if
          flsfac=lsfac
        end if

        do ir=1,nrmax
         r=(ir-1)*dr
         if(r<epslim) r=epslim
         optin(ir)=fm*(uoptcal(ir)+flsfac*uoptls(ir)-ecm)+fll1/r**2.d0
        end do

        call numerovl(ll,optin,dr,nrmax,fnumer)
        uwfsvls(ll+1,:,ijud)=fnumer(:)
       end do
      end do
      if(ilscal==0) uwfsvls(1:lmaxw+1,1:nrmax,2)=uwfsvls(1:lmaxw+1,1:nrmax,1)

      return

end subroutine wfqm
!=======================================================================



!=======================================================================
subroutine deltacalc(lmaxw,ircn,dr,hp,hm,hpd,hmd,uwfsv,smtx)
      use dims,only:ldim,nrdim
      implicit real*8(a-h,o-z)
!=======================================================================
      integer,intent(in) :: lmaxw,ircn
      real*8,intent(in) :: dr
      complex*16,intent(in) :: hp(ldim),hm(ldim),hpd(ldim),hmd(ldim),uwfsv(ldim,nrdim)

      complex*16,intent(out) :: smtx(ldim)

      complex*16 wfd

!-----------------------------------------------------------------------

      smtx(:)=(1.0d0,0.0d0)

      do l=1,lmaxw+1

       wfd=3.0d0*(uwfsv(l,ircn+1)-uwfsv(l,ircn-1))/(4.0d0*dr)  &
     &    -3.0d0*(uwfsv(l,ircn+2)-uwfsv(l,ircn-2))/(20.0d0*dr) &
     &    +1.0d0*(uwfsv(l,ircn+3)-uwfsv(l,ircn-3))/(60.0d0*dr)

       smtx(l)=(uwfsv(l,ircn)*hmd(l)-wfd*hm(l))/(uwfsv(l,ircn)*hpd(l)-wfd*hp(l))

       if(abs(smtx(l)-1.0d0)<1.0d-10) smtx(l)=(1.0d0,0.0d0)

      end do

      return

end subroutine deltacalc
!=======================================================================



!=======================================================================
subroutine wfnrm(lmaxw,nrmax,ircn,fkay,dr,hp,hm,smtxls,uwfsvls,chilud)
      use consts,only:epslim,wi
      use dims,only:ldim,nrdim
      implicit real*8 (a-h,o-z)
!=======================================================================
      integer,intent(in) :: lmaxw,nrmax,ircn
      real*8,intent(in) :: fkay,dr
      complex*16,intent(in) :: smtxls(ldim,2),uwfsvls(ldim,nrdim,2),hp(ldim),hm(ldim)
      complex*16,intent(out) :: chilud(ldim,nrdim,2)

      real*8 ra(nrdim)
      complex*16 smtxwup,smtxwdw
      complex*16 constlup(ldim),constldw(ldim)
      complex*16 pwlup,pwldw

!-----------------------------------------------------------------------

      do ir=1,nrmax
       ra(ir)=(ir-1)*dr
      end do

      do l=1,lmaxw+1
       smtxwup=smtxls(l,1)
       smtxwdw=smtxls(l,2)
       if(abs(smtxwup-1.d0)<1.d-10) smtxwup=1.0d0
       if(abs(smtxwdw-1.d0)<1.d-10) smtxwdw=1.0d0
       pwlup=(wi/2.0d0)*(hm(l)-smtxwup*hp(l))
       pwldw=(wi/2.0d0)*(hm(l)-smtxwdw*hp(l))
       constlup(l)=pwlup/uwfsvls(l,ircn,1)
       constldw(l)=pwldw/uwfsvls(l,ircn,2)
      end do

      do ir=1,nrmax
       r=(ir-1)*dr
       if(r<epslim) r=epslim
       fkr=fkay*r
       do l=1,lmaxw+1
        chilud(l,ir,1)=uwfsvls(l,ir,1)*constlup(l)
        chilud(l,ir,2)=uwfsvls(l,ir,2)*constldw(l)
       end do
      end do

      return

end subroutine wfnrm
!=======================================================================



!=======================================================================
subroutine leg_inpl(nrmax,nrgmax,lmaxw,dr,rcalg,chilud,chiludg)
      use dims,only:nrdim,ngdim,ldim
      implicit none
!=======================================================================
      integer,intent(in) :: nrmax,nrgmax,lmaxw
      real*8,intent(in) :: dr,rcalg(ngdim)
      complex*16,intent(in) :: chilud(ldim,nrdim,2)

      complex*16,intent(out) :: chiludg(ldim,ngdim,2)

      integer ir,l,i
      real*8 r,ra(nrdim),polintm

!-----------------------------------------------------------------------

      do ir=1,nrmax
       ra(ir)=(ir-1)*dr
      end do

      do i=1,2
       do ir=1,nrgmax
        r=rcalg(ir)
        do l=1,lmaxw+1
         chiludg(l,ir,i)=                                                 &
     &    dcmplx(polintm(ra(1:nrmax),dble(chilud(l,1:nrmax,i)),nrmax,r,3) &
     &          ,polintm(ra(1:nrmax),imag(chilud(l,1:nrmax,i)),nrmax,r,3))
        end do
       end do
      end do

      return

end subroutine leg_inpl
!=======================================================================



!=======================================================================
subroutine makechil(lmaxw,nrgmax,fkay,sigma,rcalg,chiludg,childkruug,childkrudg)
      use consts,only:epslim,wi
      use dims,only:ldim,ngdim
      implicit none
!=======================================================================
      integer,intent(in) :: lmaxw,nrgmax
      real*8,intent(in) :: fkay,sigma(ldim),rcalg(ngdim)
      complex*16,intent(in) :: chiludg(ldim,ngdim,2)
      complex*16,intent(out) :: childkruug(ldim,ngdim),childkrudg(ldim,ngdim)

      integer ir,l
      real*8 r,fkr
      complex*16 factlsv0(ldim),factlsv1(ldim)

!-----------------------------------------------------------------------

      do l=1,lmaxw+1
       factlsv0(l)=wi**(l-1)*exp(wi*sigma(l))
       factlsv1(l)=wi**(l-1)*sqrt((l-1.d0)*l)*exp(wi*sigma(l))
      end do

      do ir=1,nrgmax
       r=rcalg(ir)
       fkr=fkay*r
       do l=1,lmaxw+1
        childkruug(l,ir)=factlsv0(l)*(l*chiludg(l,ir,1)+(l-1)*chiludg(l,ir,2))
        childkrudg(l,ir)=factlsv1(l)*(chiludg(l,ir,1)-chiludg(l,ir,2))
       end do
       childkruug(1:lmaxw+1,ir)=childkruug(1:lmaxw+1,ir)/fkr
       childkrudg(1:lmaxw+1,ir)=childkrudg(1:lmaxw+1,ir)/fkr
      end do

end subroutine makechil
!=======================================================================



!=======================================================================
subroutine make_pfacsvg(nrmax,nrgmax,beta,dr,fmu,uopt,fnl,rcalg,pfacsvg)
      use consts,only:ac,hc,epslim
      use dims,only:ngdim,nrdim
      implicit none
!=======================================================================
      integer,intent(in) :: nrmax,nrgmax
      real*8,intent(in) :: beta,dr,fmu,rcalg(ngdim)
      complex*16,intent(in) :: uopt(nrdim),fnl(nrdim)

      complex*16,intent(out) :: pfacsvg(ngdim)

      integer ir
      real*8 r,ra(nrdim),facp,polintm
      complex*16 pfacsv(nrdim)

!-----------------------------------------------------------------------

      do ir=1,nrmax
       ra(ir)=(ir-1)*dr
      end do

      if(beta>epslim) then
        facp=fmu*ac*beta**2.0d0/(2.0d0*hc**2.0d0)
        do ir=1,nrmax
         pfacsv(ir)=(1.0d0-facp*uopt(ir))**(-0.5d0)
  !---another choice of the perey correction used in
  !   timofeyuk and johnson, phys. rev. c 87, 064610 (2013).----
  !      pfacsv(ir)
  !        =exp(fmu*ac*beta**2.0d0*uopt(ir)/(4.0d0*hc**2.0d0))
        end do
       else if(beta<-epslim) then
        do ir=1,nrmax
         pfacsv(ir)=fnl(ir)
        end do
      end if

      do ir=1,nrgmax
       r=rcalg(ir)
       pfacsvg(ir)=dcmplx(polintm(ra(1:nrmax),dble(pfacsv(1:nrmax)),nrmax,r,3), &
     &                    polintm(ra(1:nrmax),imag(pfacsv(1:nrmax)),nrmax,r,3))
        end do

      return

end subroutine make_pfacsvg
!=======================================================================



!=======================================================================
subroutine wfout(lminw,lmaxw,nthgmax,nphgmax,nrgmax,nzplmax,ilscal &
     &          ,fkay,thtk,phik,pm,thcalg,phcalg,rcalg,pltbl,dzpl  &
     &          ,wfls,childkruug,childkrudg)
      use consts,only:epslim,pi,wi
      use dims,only:ldim,ngdim,nrdim,nzpldim
      implicit real*8 (a-h,o-z)
!=======================================================================
      integer,intent(in) :: lminw,lmaxw,nthgmax,nphgmax &
     &                     ,nrgmax,nzplmax,ilscal
      real*8,intent(in) :: fkay,thtk,phik,pm,thcalg(ngdim),phcalg(ngdim),rcalg(ngdim) &
     &                    ,pltbl(ldim,nzpldim),dzpl

      complex*16,intent(out) :: wfls(2,2,ngdim,ngdim,ngdim)
      complex*16,intent(inout) :: childkruug(ldim,ngdim),childkrudg(ldim,ngdim)

      integer iargsv(ngdim,ngdim)
      real*8 plcalsv(ngdim,ngdim,ldim)
      real*8 phiftbl(ngdim,ngdim)
      real*8 dfunctbl(nthgmax,nphgmax,lmaxw+1)
      real*8 yl1tbl(nthgmax,nphgmax,lmaxw+1)
      real*8 cosgtbl(nthgmax,nphgmax)
      real*8 ylmtbl(nphgmax,lmaxw+1,lmaxw+1)
      complex*16 wf(ngdim,ngdim,ngdim)
      complex*16 wfl
      complex*16 wfluu,wflud
      complex*16 wfluusum,wfludsum

!-----------------------------------------------------------------------

      wfls(:,:,1:nrgmax,1:nthgmax,1:nphgmax)=(0.d0,0.d0)

      if(ilscal==0) then

        wf(1:nrgmax,1:nthgmax,1:nphgmax)=(0.d0,0.d0)

        if(pm<0.d0) then
          do l=lminw+1,lmaxw+1
           childkruug(l,1:nrgmax)=childkruug(l,1:nrgmax)*(-1.d0)**(l-1)
          end do
        end if

        do ithg=1,nthgmax
         tht=thcalg(ithg)
         do iphg=1,nphgmax
          phi=phcalg(iphg)
          call rotate(tht,phi,thtk,phik,thtf,phif)
          cosg=cos(thtf)
          iarg=int((cosg+1.0d0)/dzpl+1.01)
          iargsv(ithg,iphg)=iarg
          cosg0=(iarg-1)*dzpl-1.0d0
          fac=(cosg-cosg0)/dzpl
          do l=lminw+1,lmaxw+1
           if(iarg==nzplmax) then
             plcalsv(ithg,iphg,l)=pltbl(l,iarg)
            else
             plcalsv(ithg,iphg,l)=(pltbl(l,iarg+1)-pltbl(l,iarg))*fac+pltbl(l,iarg)
           end if
          end do
         end do
        end do

        do ir=1,nrgmax
         r=rcalg(ir)
         if(r==0.0d0) r=epslim
         fkr=fkay*r
         do ithg=1,nthgmax
          tht=thcalg(ithg)
          do iphg=1,nphgmax
           phi=phcalg(iphg)

           iarg=iargsv(ithg,iphg)

           do l=lminw+1,lmaxw+1
            plcal=plcalsv(ithg,iphg,l)
            wfl=childkruug(l,ir)*plcal
            wf(ir,ithg,iphg)=wf(ir,ithg,iphg)+wfl
           end do

          end do
         end do
        end do

        wfls(1,1,1:nrgmax,1:nthgmax,1:nphgmax)=wf(1:nrgmax,1:nthgmax,1:nphgmax)

      else if(ilscal==1) then

       if(pm<0.d0) then ! conversion from outgoing to incoming boundary condition
         do l=lminw+1,lmaxw+1
          childkruug(l,1:nrgmax)=childkruug(l,1:nrgmax)*(-1.d0)**(l-1)
          childkrudg(l,1:nrgmax)=childkrudg(l,1:nrgmax)*(-1.d0)**(l-1)
         end do
       end if

       phiftbl(:,:)=0.d0
       do ithg=1,nthgmax
        tht=thcalg(ithg)
        do iphg=1,nphgmax
         phi=phcalg(iphg)
         call rotate(tht,phi,thtk,phik,thtf,phif)
         phiftbl(ithg,iphg)=phif
         cosg=cos(thtf)
         cosgtbl(ithg,iphg)=cosg
         iarg=int((cosg+1.0d0)/dzpl+1.01)
         iargsv(ithg,iphg)=iarg
         cosg0=(iarg-1)*dzpl-1.0d0
         fac=(cosg-cosg0)/dzpl
         do l=lminw+1,lmaxw+1
          if(iarg==nzplmax) then
            plcalsv(ithg,iphg,l)=pltbl(l,iarg)
           else
            plcalsv(ithg,iphg,l)=(pltbl(l,iarg+1)-pltbl(l,iarg))*fac+pltbl(l,iarg)
          end if
         end do
        end do
       end do

       yl1tbl(:,:,:)=0.d0
       dfunctbl(:,:,:)=0.d0
       if(lmaxw>0) then
         do ithg=1,nthgmax
          call ylmph0(lmaxw+1,lmaxw+1,nphgmax,lmaxw,1,1,nphgmax,cosgtbl(ithg,:),ylmtbl)
          yl1tbl(ithg,1:nphgmax,1:lmaxw+1)=ylmtbl(1:nphgmax,1:lmaxw+1,2)
         end do

         dfunctbl(1:nthgmax,1:nphgmax,1)=1.d0
         do l=lminw+2,lmaxw+1
          dfunctbl(1:nthgmax,1:nphgmax,l) &
           =sqrt(4.d0*pi/(2.d0*(l-1.d0)+1.d0))*yl1tbl(1:nthgmax,1:nphgmax,l)
         end do
       end if

       do ir=1,nrgmax
        r=rcalg(ir)
        if(r==0.0d0) r=epslim
        fkr=fkay*r
        do ithg=1,nthgmax
         tht=thcalg(ithg)
         do iphg=1,nphgmax
          phi=phcalg(iphg)

          iarg=iargsv(ithg,iphg)

          wfluusum=(0.d0,0.d0)
          wfludsum=(0.d0,0.d0)

          do l=lminw+1,lmaxw+1
           plcal=plcalsv(ithg,iphg,l)
           dfunccal=dfunctbl(ithg,iphg,l)
           wfluu=childkruug(l,ir)*plcal
           wflud=childkrudg(l,ir)*dfunccal
           wfluusum=wfluusum+wfluu
           wfludsum=wfludsum+wflud
          end do

          phif=phiftbl(ithg,iphg)
          wfls(1,1,ir,ithg,iphg)=wfluusum
          wfls(1,2,ir,ithg,iphg)=wfludsum*exp(pm*wi*phif)
          wfls(2,1,ir,ithg,iphg)=-wfludsum*exp(-pm*wi*phif)
          wfls(2,2,ir,ithg,iphg)=wfluusum

         end do
        end do
       end do

      end if

      return

end subroutine wfout
!=======================================================================



!=======================================================================
subroutine tmtx(dr,fkayi,thtkg,amasspi,aa,ylm0,wfmtp               &
     &         ,ffrg,nrmax,ls,thcalg,thweig,nthgmax,phweig,nphgmax &
     &         ,rcalg,rweig,nrgmax,iex,io,itmdcal,allsumm,mvalmin,ylmphi2tbl)
      use consts
      use dims
      use kibmod
      implicit real*8(a-h,o-z)
!=======================================================================
      real*8 fkayi(3),thtkg(3),amasspi(3)
      real*8 ylm0(nthdim,ldim,mdim)
      real*8 ffrg(ngdim)
      complex*16 wfmtp(ngdim,ngdim,ngdim)
      complex*16 thtintg(ngdim)
      complex*16 rintgrnd(ngdim)
      complex*16 deltam(-mdim:mdim,ngdim),delta(ngdim)
      complex*16 allsumm(-mdim:mdim)
      complex*16 allsum,gaus,gaustht,funcinteg1,faccor
      real*8 thcalg(ngdim),thweig(ngdim)
      real*8 phweig(ngdim)
      real*8 rcalg(ngdim),rweig(ngdim)
      complex*16 ylmphi2tbl(-ldim:ldim,ngdim)

!-----------------------------------------------------------------------

      allmsum=0.0d0

      deltam(:,:)=(0.0d0,0.0d0)

      amassc=amasspi(3)
      if(iex==1) amassc=amasspi(2)

      do mval=mvalmin,ls

       thtintg=(0.0d0,0.0d0)
       allsum=(0.0d0,0.0d0)

       do ir=1,nrgmax
        r=rcalg(ir)
        fkdr=fkayi(1)*cos(thtkg(1))*r

        gaus=(0.0d0,0.d0)
        sigr=1.0d0

        do ithg=1,nthgmax
         gaustht=(0.d0,0.d0)
         do iphg=1,nphgmax
          funcinteg1=wfmtp(ir,ithg,iphg)*ylmphi2tbl(mval,iphg)
          gaustht=gaustht+funcinteg1*phweig(iphg)
         end do
         fkdrcos=fkdr*cos(thcalg(ithg))
         faccor=exp(-wi*fkdrcos*amassc/aa)
         gaus=gaus+gaustht*thweig(ithg)*ylm0(ithg,ls+1,abs(mval)+1) &
     &        *sin(thcalg(ithg))*faccor
        end do

        thtintg(ir)=gaus*sigr

        rintgrnd(ir)=thtintg(ir)*ffrg(ir)*r**2.0d0
        allsum=allsum+rintgrnd(ir)*rweig(ir)
       end do

       allsumm(mval)=allsum

       if(itmdcal==1) then
         deltam(mval,1:nrgmax)=conjg(allsum)*rintgrnd(1:nrgmax)
       end if
       allsum2=abs(allsum)**2.0d0

       if(mval/=0 .and. mvalmin==0) allsum2=2.0d0*allsum2

       allmsum=allmsum+allsum2

      end do

      if(itmdcal==1) &
     &  call tdenssum(mvalmin,ls,nrmax,dr,deltam,allmsum,delta,nrgmax,rcalg,rweig,io)

      return

end subroutine tmtx
!=======================================================================



!=======================================================================
subroutine kinemaelm(iang,fk1,e1,th1,ph1,fk2,e2                  &
     &              ,th2,ph2,aa,fk0,th0,ph0,fm0,fmsp,fm1,fm2,iex &
     &              ,fkapf,fkapi,thnn,fkapfe,fkapie,thnne        &
     &              ,thkapie,thkapfe,phkapie,phkapfe             &
     &              ,ennin,beteg,fmol2,fkaponsh,fmurel)
      use consts
      use nntbl,only:ionsh
      implicit real*8 (a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

!     fkapi and fkapf is initial and final relative momenta of nn system
!     in G-frame. fkapie and fkapfe are those in c.m. frame in nn system.

      call vecxyz(iang,fk0,th0,ph0, fk0x,fk0y,fk0z)
      call vecxyz(iang,fk1,th1,ph1, fk1x,fk1y,fk1z)
      call vecxyz(iang,fk2,th2,ph2, fk2x,fk2y,fk2z)

      coefa=1.0d0
      coefi1=fmsp/(fm0+fmsp)
      coefi2=fm0/(fm0+fmsp)
      coeff1=fm2/(fm1+fm2)
      coeff2=fm1/(fm1+fm2)
      fkapfx=coeff1*fk1x-coeff2*fk2x
      fkapfy=coeff1*fk1y-coeff2*fk2y
      fkapfz=coeff1*fk1z-coeff2*fk2z

      fk0efx=coefa*fk0x
      fk0efy=coefa*fk0y
      fk0efz=coefa*fk0z
      fkspx=fk1x+fk2x-fk0efx
      fkspy=fk1y+fk2y-fk0efy
      fkspz=fk1z+fk2z-fk0efz

      fkapix=coefi1*fk0efx-coefi2*fkspx
      fkapiy=coefi1*fk0efy-coefi2*fkspy
      fkapiz=coefi1*fk0efz-coefi2*fkspz

      fkapf=sqrt(fkapfx**2+fkapfy**2+fkapfz**2)
      fkapi=sqrt(fkapix**2+fkapiy**2+fkapiz**2)

      sprokk=fkapfx*fkapix+fkapfy*fkapiy+fkapfz*fkapiz
      fkk=fkapf*fkapi
      if(fkk<epslim) fkk=epslim
      costhnn=sprokk/fkk
      thnn=acos2(costhnn)/piad

!-----G frame to the cm frame of the colliding two particles

      fac=hc/(e1+e2)
      betegx=fac*(fk1x+fk2x)
      betegy=fac*(fk1y+fk2y)
      betegz=fac*(fk1z+fk2z)
      beteg=sqrt(betegx**2+betegy**2+betegz**2)
      gameg=1.0d0/sqrt(1.0d0-beteg**2)

      call ltrgen(fk1x,fk1y,fk1z,e1,betegx,betegy,betegz,beteg,gameg &
     &           ,fk1ex,fk1ey,fk1ez,e1e)
      call ltrgen(fk2x,fk2y,fk2z,e2,betegx,betegy,betegz,beteg,gameg &
     &           ,fk2ex,fk2ey,fk2ez,e2e)

      e0ef=sqrt((hc*coefa*fk0)**2+fm0**2)
      fksp=sqrt(fkspx**2+fkspy**2+fkspz**2)
      esp =sqrt((hc*fksp)**2+fmsp**2)

      call ltrgen(fk0efx,fk0efy,fk0efz,e0ef,betegx,betegy,betegz &
     &           ,beteg,gameg, fk0ex,fk0ey,fk0ez,e0e)
      call ltrgen(fkspx,fkspy,fkspz,esp,betegx,betegy,betegz &
     &           ,beteg,gameg, fkspex,fkspey,fkspez,espe)

      fkapfex=coeff1*fk1ex-coeff2*fk2ex
      fkapfey=coeff1*fk1ey-coeff2*fk2ey
      fkapfez=coeff1*fk1ez-coeff2*fk2ez
      fkapiex=coefi1*fk0ex-coefi2*fkspex
      fkapiey=coefi1*fk0ey-coefi2*fkspey
      fkapiez=coefi1*fk0ez-coefi2*fkspez

      fkapfe=sqrt(fkapfex**2+fkapfey**2+fkapfez**2)
      fkapie=sqrt(fkapiex**2+fkapiey**2+fkapiez**2)
      sprokke=fkapfex*fkapiex+fkapfey*fkapiey+fkapfez*fkapiez
      fkke=fkapfe*fkapie
      if(fkke<epslim) fkke=epslim
      costhnne=sprokke/fkke
      thnne=acos2(costhnne)/piad
      if(fkapie<epslim) fkapie=epslim
      if(fkapfe<epslim) fkapfe=epslim
      thkapie=acos2(fkapiez/fkapie)/piad
      thkapfe=acos2(fkapfez/fkapfe)/piad

      phkapie=sign(1.d0,fkapiey)*acos2(fkapiex/sqrt(fkapiex**2.d0+fkapiey**2.d0))
      if(fkapiex==0.d0 .and. fkapiey==0.d0) phkapie=0.d0
      if(phkapie<0.d0) phkapie=phkapie+2.d0*pi
      phkapie=phkapie/piad

      phkapfe=sign(1.d0,fkapfey)*acos2(fkapfex/sqrt(fkapfex**2.d0+fkapfey**2.d0))
      if(fkapfex==0.d0 .and. fkapfey==0.d0) phkapfe=0.d0
      if(phkapfe<0.d0) phkapfe=phkapfe+2.d0*pi
      phkapfe=phkapfe/piad

      fmu=fm1*fm2/(fm1+fm2)

      fmol2=e1e*e2e*e0e*espe/(e1*e2*e0ef*esp)

      fmureli=e0e*esp/(e0e+esp)
      fmurelf=e1e*e2e/(e1e+e2e)

      if(iex==0) then
        bet2x=fk2x*hc/e2
        bet2y=fk2y*hc/e2
        bet2z=fk2z*hc/e2
        bet2=sqrt(bet2x**2+bet2y**2+bet2z**2)
        gam2=1.0d0/sqrt(1.0d0-bet2**2)

        call ltrgen(fk1x,fk1y,fk1z,e1,bet2x,bet2y,bet2z,bet2,gam2,fk12x,fk12y,fk12z,e12)
        enninf=e12-fm0

       else if(iex==1) then

        bet1x=fk1x*hc/e1
        bet1y=fk1y*hc/e1
        bet1z=fk1z*hc/e1
        bet1=sqrt(bet1x**2+bet1y**2+bet1z**2)
        gam1=1.0d0/sqrt(1.0d0-bet1**2)

        call ltrgen(fk2x,fk2y,fk2z,e2,bet1x,bet1y,bet1z,bet1,gam1,fk21x,fk21y,fk21z,e21)
        enninf=e21-fm0

      end if

      bet2x=fkspx*hc/esp
      bet2y=fkspy*hc/esp
      bet2z=fkspz*hc/esp
      bet2=sqrt(bet2x**2+bet2y**2+bet2z**2)
      gam2=1.0d0/sqrt(1.0d0-bet2**2)

      call ltrgen(fk0efx,fk0efy,fk0efz,e0ef,bet2x,bet2y,bet2z,bet2,gam2 &
     &           ,fk0spx,fk0spy,fk0spz,e0sp)
      ennini=e0sp-fm0

!     setting E_{NN} and kappa_{NN} in the on-shell approx.

!     BE CAREFUL!
!     E_{NN} is defined in L-frame of NN 2-body system, while kappa_{NN} is in t-frame.

      if(ionsh==1) then ! final state prescription
        ennin=enninf/(fm0/ac)
        fkaponsh=fkapfe
        fmurel=fmurelf
       else if(ionsh==2) then ! initial state prescription
        ennin=ennini/(fm0/ac)
        fkaponsh=fkapie
        fmurel=fmureli
       else if(ionsh==3) then ! energy average prescription
        ennin=(enninf+ennini)/2.d0/(fm0/ac)
        fmurel=(fmureli+fmurelf)/2.d0
        fkaponsh=sqrt((fkapie**2.d0+fkapfe**2.d0)/2.d0)
       else if(ionsh==4) then ! momentum average prescription
        fki=sqrt(2.d0*fmureli*ac*ennini)/hc
        fkf=sqrt(2.d0*fmurelf*ac*enninf)/hc
        fkav=(fki+fkf)/2.d0
        fmurel=(fmureli+fmurelf)/2.d0
        ennin=(hc*fkav)**2.d0/(2.d0*fmurel*ac)/(fm0/ac)
        fkaponsh=(fkapfe+fkapie)/2.d0
      end if

      return

end subroutine kinemaelm
!=======================================================================



!=======================================================================
subroutine sigmal(eta,lmax1,sigma,exsig)
      use dims
      implicit real*8(a-h,o-z)
!=======================================================================
      complex*16 exsig(ldim)
      real*8 sigma(ldim)

!-----------------------------------------------------------------------

      lmax0=lmax1-1
      if(lmax1>ldim) go to 9000
      if(eta>=10.0)  go to 10
      eta2=eta*eta
      eta2a=eta+eta
      eta6=eta2+16.0d0
      sigma0=-(eta/(12.0*eta6))*(1.0+(eta2-48.0)/(30.0*(eta6**2)) &
     &        +((eta2-160.0)*eta2+1280.0)/(105.0*(eta6**4)))      &
     &        -eta+(eta/2.0)*dlog(eta6)+3.5*datan(0.25*eta)       &
     &        -(datan(eta)+datan(0.5*eta)+datan(eta/3.0))
      go to 11
   10 einv1=1.0/eta
      einv2=einv1*einv1
      einv3=einv1*einv2
      einv5=einv3*einv2
      einv7=einv5*einv2
      einv9=einv7*einv2
      sigma0=0.7853981634d0+eta*dlog(eta)-eta              &
     &       -(0.08333333333d0*einv1+0.00277777777d0*einv3 &
     &       +0.00079365079d0*einv5+0.00059523810d0*einv7  &
     &       +0.00084175084d0*einv9)
   11 modtpi=nint(sigma0/6.2831853071796d0)
      sigmaz=sigma0-6.2831853071796d0*modtpi
      sigma(1)=sigma0
      exsig(1)=dcmplx(dcos(sigmaz),dsin(sigmaz))
      if(lmax0<=0) return
      do 1 ll=1,lmax0
      ll1=ll+1
      sigma(ll1)=sigma(ll)+datan(eta/ll)
      modtpi=nint(sigma(ll1)/6.2831853071796d0)
      sigmaz=sigma(ll1)-6.2831853071796d0*modtpi
      exsig(ll1)=dcmplx(dcos(sigmaz),dsin(sigmaz))
    1 continue
      return
 9000 continue
      write(10,9100) lmax1,ldim
 9100 format(1x/10x,'*** error *** dimension over in sigmal'/ &
     &   20x,' : lmax1(=',i4,' > ldim(=',i4,') '/)
      stop 1
      return

end subroutine sigmal
!=======================================================================



!=======================================================================
function  plfunc(l,z)
      implicit real*8(a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      plfunc=0.d0
      if(l>200) go to 50
      if(l>10) go to 21
      ll=l+1
      go to (10,11,12,13,14,15,16,17,18,19,20),ll
   10 continue
      plfunc=1.0
      return
   11 continue
      plfunc=z
      return
   12 continue
      plfunc=(3.0*z**2 - 1.0)/2.0
      return
   13 continue
      plfunc=(5*z**2 - 3  )*z/2.0
      return
   14 continue
      z2=z**2
      plfunc=(35*z2**2 - 30*z2 + 3)/8.0
      return
   15 continue
      z2=z**2
      plfunc=(63*z2**2 - 70*z2 + 15)*z/8.0
      return
   16 continue
      z2=z**2
      plfunc=(231*z2**3 - 315*z2**2 +105*z2 -5)/16.0
      return
   17 continue
      z2=z**2
      plfunc=(429*z2**3 - 693*z2**2 + 315*z2 -35)*z/16.0
      return
   18 continue
      z2=z**2
      plfunc=(6435*z2**4 - 12012*z2**3 + 6930*z2**2 - 1260*z2 +35)/128.0
      return
   19 continue
      z2=z**2
      plfunc=(12155*z2**4 - 25740*z2**3 + 18018*z2**2 - 4620*z2 + 315)*z/128.0
      return
   20 continue
      z2=z**2
      plfunc=(230945*z2**5 - 546975*z2**4 + 450450*z2**3 - 150150*z2**2 &
     &        + 17325*z2 - 315)/1280.0
      return
   21 continue
      z2=z**2
      plpre2=(12155*z2**4 - 25740*z2**3 + 18018*z2**2 - 4620*z2 + 315)*z/128.0
      plpre1=(230945*z2**5 - 546975*z2**4 + 450450*z2**3 - 150150*z2**2 &
     &        + 17325*z2 - 315)/1280.0
      ilmax=l-10
      lcal=11
      do il=1,ilmax
       plfunc=( (2.0d0*lcal-1.0d0)*z*plpre1 - (lcal-1.0d0)*plpre2 )/(lcal*1.0d0)
       plpre2=plpre1
       plpre1=plfunc
       lcal=lcal+1
      end do
      return

   50 continue
      write(10,900) l
  900 format(/1x,'*******  error in plfunc --- no definition for l>200',': l=',i5)
      stop 1

end function plfunc
!=======================================================================



!=======================================================================
!---- derivative of Legendre polynomial p_l(z)
function  dplfunc(l,z)
      implicit real*8(a-h,o-z)
!=======================================================================

!-----------------------------------------------------------------------

      if(l==0) then
        dplfunc=0.0d0
        return
      end if

      if(abs(z)==1.0d0) then
        dplfunc=z**(l+1)*l*(l+1)/2.0d0
       else
        dplfunc=l*(z*plfunc(l,z)-plfunc(l-1,z))/(z**2.0d0-1.0d0)
      end if

      return

end function dplfunc
!=======================================================================



!=======================================================================
!******************************************
!********    ylmphi2    *******************
!******************************************
function ylmphi2(m,phi)
      implicit real*8(a-h,o-z)
!=======================================================================
      complex*16 ylmphi2,wi

!-----------------------------------------------------------------------

      wi=(0.0d0,1.0d0)

      cosphi=cos(phi)
      sinphi=sin(phi)
      ylmphi2=(cosphi+wi*sinphi)**m
      if(m<0) ylmphi2=(-1.0d0)**(-m)*ylmphi2

      return

end function ylmphi2
!=======================================================================



!=======================================================================
!**
!**        ***** clockm *****
!**          function etime ......... cpu timer for epcf77 (sun)
!**          subroutine cpu_time .... cpu timer for fortran95
!**
subroutine clockm(itime)
!  -- etime
!      dimension time(2)
!      t=etime(time)
!      t=etime_(time)
!      itime=nint(time(1)*1000)
!  -- cpu_time
      call cpu_time(ctime)
      itime=nint(ctime*1000)
      return

end subroutine clockm
!=======================================================================



!=======================================================================
subroutine tdenssum(mvalmin,ls,nrmax,dr,deltam,allmsum,delta,nrgmax,rcalg,rweig,io)
      use consts
      use dims
      use kibmod
      implicit real*8(a-h,o-z)
!=======================================================================
      complex*16 checkmsum,deltainpol
      complex*16 deltam(-mdim:mdim,ngdim),delta(ngdim)
      real*8     rcalg(ngdim),rweig(ngdim)

!-----------------------------------------------------------------------

      delta(:)=(0.0d0,0.0d0)
      checkmsum=(0.0d0,0.0d0)
      fmdens=0.0d0

      write(kibtmd,'(/,a35)') 'transition matrix density output'
      write(kibtmd,'(a8,i10)') 'config#:',io
      write(kibtmd,'(3a15)') 'r','dble(delta)','imag(delta)'
      do ir=1,nrgmax
       r=rcalg(ir)
       if(mvalmin==0) then
         delta(ir)=deltam(0,ir)+2.0d0*sum(deltam(1:ls,ir))
        else if(mvalmin==-ls) then
         delta(ir)=sum(deltam(-ls:ls,ir))
       end if

      end do

      do ir=1,nrmax
       r=(ir-1)*dr
       deltainpol=                                                         &
     &   dcmplx(polintm(rcalg(1:nrgmax),dble(delta(1:nrgmax)),nrgmax,r,3), &
     &          polintm(rcalg(1:nrgmax),imag(delta(1:nrgmax)),nrgmax,r,3) )
       write(kibtmd,'(3e15.5)') r,deltainpol
      end do

      do ir=1,nrgmax
       checkmsum=checkmsum+delta(ir)*rweig(ir)
      end do

      diff=dble(checkmsum)/allmsum
      if(abs(1.0d0-diff)>chkcons) then
        write(kibout,*) 'delta_r error in tdenssum'
        write(kibout,'(3a15)') 'checkmsum','allmsum','error'
        write(kibout,'(3e15.5)') dble(checkmsum),allmsum,abs(1.0d0-diff)
        stop 1
       else if(abs(aimag(checkmsum)/dble(checkmsum))>1.0d-3) then
        write(kibout,*) 'delta_i error in tdenssum'
        write(kibout,'(3a15)') 'real_checksum','imag_checksum','error'
        write(kibout,'(3e15.5)') dble(checkmsum),aimag(checkmsum), &
     &                           aimag(checkmsum)/dble(checkmsum)
        stop 1
      end if

      return

end subroutine tdenssum
!=======================================================================



!=======================================================================
subroutine rotate(thtri,phiri,thtki,phiki,thtrf,phirf)
      use consts,only:pi,epslim
!=======================================================================
!
!     function to rotate the coordinate.
!     this is related to euler angles with
!     active rotation, z-y-z convention.
!
!     for any r and k,
!     wfi(r,thtri,phiri,k,thtki,phiki)=wff(r,thtrf,phirf,k,thtkf,phikf).
!
!     all angles should be given in radian.
!
!     for practical use, input the direction of r and k into
!     thtri,phiri,thtki and phiki.
!     the subrutine reterns thtrf and phirf, which are the rotation
!     angles from k to r, i.e., coordinate in helicity frame.
!
!     input:
!      thtki, phiki: direction of the momentum of the scatttering wave.
!      thtri, phiri: the angles of the coordinate r that shows the
!                    position of the wave function you want.
!     output:
!      thtrf, phirf: the rotation angles from the direction of k to r,
!                    i.e., the coordinate in helicity frame.
!
!=======================================================================
      implicit real*8 (a-h,o-z)
!=======================================================================

      xi=sin(thtri)*cos(phiri)
      yi=sin(thtri)*sin(phiri)
      zi=cos(thtri)

      ct=cos(thtki)
      st=sin(thtki)
      cp=cos(phiki)
      sp=sin(phiki)

      xf=ct*cp*xi+ct*sp*yi-st*zi
      yf=  -sp*xi+cp*yi
      zf=cp*st*xi+st*sp*yi+ct*zi

      if(abs(xf)<epslim) xf=epslim

      thtrf=0.0d0
      phirf=0.0d0

      if(abs(xf)<=epslim .and. abs(yf)<=epslim .and. abs(zf)<=epslim) then
        thtrf=0.d0
       else
        thtrf=acos2(zf)
      end if
      if(abs(xf)<=epslim .and. abs(yf)<=epslim) then
        phirf=0.d0
       else
        phirf=sign(1.d0,yf)*acos2(xf/sqrt(xf**2.d0+yf**2.d0))
      end if
      if(phirf<0.0d0) phirf=phirf+2.d0*pi

      return

end subroutine rotate
!=======================================================================



!=======================================================================
subroutine onshnnread(nthqdim,nphqdim,nq0dim,nthq0dim,nphq0dim,gma &
     &               ,nq0mn,nq0mx,nth0mn,nth0mx,nph0mn,nph0mx,nthmn,nthmx,nphmn,nphmx)
      use kibmod,only:kibelm
      implicit real*8 (a-h,o-z)
!=======================================================================
      complex*16 gma(nthqdim,nphqdim,nq0dim,nthq0dim,nphq0dim,16)

!-----------------------------------------------------------------------

      iline=0
      do ! for file length check
       do iq0=nq0mn,nq0mx
        do ith0=nth0mn,nth0mx
         do iph0=nph0mn,nph0mx
          do ith=nthmn,nthmx
           do iph=nphmn,nphmx
            do ispin=1,16
             read(kibelm,*,end=999) ampr,ampi
             iline=iline+1
             gma(ith,iph,iq0,ith0,iph0,ispin)=dcmplx(ampr,ampi)
            end do ! ispin
           end do ! iph
          end do ! ith
         end do ! iph0
        end do ! ith0
       end do ! iq0
      end do

 999  continue
      ilinetot=(nq0mx-nq0mn+1)*(nth0mx-nth0mn+1)*(nph0mx-nph0mn+1) &
     &        *(nthmx-nthmn+1)*(nphmx-nphmn+1)*16
      if(ilinetot/=iline) then
        write(*,'(a,i0,a)') 'ERROR: total line numbers of input # ',kibelm, &
     &                      ' is inconsistent with header information'
        write(*,'(a,i0)') 'lines expected = ',ilinetot+3
        write(*,'(a,i0)') 'actual lines read-in = ',iline+3
        stop 1
      end if

      return

end subroutine  onshnnread
!=======================================================================



!=======================================================================
!**
!**       ******** clebsh ********
!**
function clebsh(fj1,fm1,fj2,fm2,fj3,fm3)
      use angmom,only:faclog,memoam
      implicit  real*8(a-h,o-z)
!**
!**   a function for clebsch-gordan coefficients with arbitrary arguments
!**   an original fortran program code of function type with integer
!**   arguments by t.tamura is revised as function type with floating
!**   point arguments by m.sakakura
!**   input is (j1,m1,j2,m2/j,m) floating point arguments
!**
      if(memoam >= nint(fj1+fj2+fj3)+2) go to 2
      memoam=nint(fj1+fj2+fj3)+2
      faclog(1)=0.0d0
      faclog(2)=0.0d0
      fn=1.0d0
      do 10 i=3,nint(fj1+fj2+fj3)+2
      fn=fn+1.0d0
   10 faclog(i)=faclog(i-1)+dlog(fn)
    2 clebsh=0.0d0
      if(fj1) 506,30,30
   30 if(fj2) 506,31,31
   31 if(fj3) 506,32,32
   32 ia = nint(2.0*fj1+0.01)
      ib = nint(2.0*fj2+0.01)
      ic = nint(2.0*fj3+0.01)
      if (fm1) 41,42,43
   41 id = nint(2.0*fm1 -0.01)
      go to 50
   42 id = 0
      go to 50
   43 id=nint(2.0*fm1+ 0.01)
   50 if (fm2) 51,52,53
   51 ie = nint(2.0*fm2- 0.01)
      go to 60
   52 ie = 0
      go to 60
   53 ie = nint(2.0*fm2 + 0.01)
   60 if (fm3) 61,62,63
   61 if = nint(2.0*fm3 - 0.01)
      go to 70
   62 if = 0
      go to 70
   63 if = nint(2.0*fm3 + 0.01)
   70 if(id+ie-if) 500,105,500
  105 k1=ia+ib+ic
      if(k1-2*(k1/2)) 501,110,501
  110 k1=ia+ib-ic
      k2=ic-iabs (ib-ia)
      k3= min0 (k1,k2)
      if(k3) 502,130,130
  130 if((-1)**(ib+ie)) 503,503,140
  140 if((-1)**(ic+if)) 503,503,150
  150 if(ia-iabs (id)) 505,152,152
  152 if(ib-iabs (ie)) 505,154,154
  154 if(ic-iabs (if)) 505,160,160
  160 if(ia) 506,175,165
  165 if(ib) 506,175,170
  170 if(ic) 506,180,250
  175 clebsh=1.0
      return
  180 fb=ib+1
      clebsh=((-1.0)**((ia-id)/2))/dsqrt(fb)
      go to 1000
  250 fc2=ic+1
      iabcp=(ia+ib+ic)/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      sqfclg=0.5*(dlog(fc2)-faclog(iabcp+1)          &
     &       +faclog(iabc)+faclog(icab)+faclog(ibca) &
     &       +faclog(iapd)+faclog(iamd)+faclog(ibpe) &
     &       +faclog(ibme)+faclog(icpf)+faclog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi= max0 (0,nzmic2,nzmic3)+1
      nzmx= min0 (iabc,iamd,ibpe)
      s1=(-1.0)**(nzmi-1)
      if(nzmx<nzmi) go to 1000
      do 400 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2) &
     &             -faclog(nzt3)-faclog(nzt4)-faclog(nzt5)
      ssterm=s1*dexp(termlg)
      clebsh=clebsh+ssterm
      s1=-s1
  400 continue
      return
  500 continue
!-----when various conditions for a cg-coefficient does not hold,
!-----this routine makes clebsh=0 (see above).
!!    write(kibout,600) fm1,fm2,fm3
      return
  501 continue
!!    write(kibout,601) fj1,fj2,fj3
      return
  502 continue
!!    write(kibout,602) fj1,fj2,fj3
      return
  503 continue
!!    write(kibout,603) fj1,fm1,fj2,fm2,fj3,fm3
      return
  505 continue
!!    write(kibout,605) fj1,fm1,fj2,fm2,fj3,fm3
      return
  506 continue
!!    write(kibout,606) fj1,fj2,fj3
      return
!  600 format(' ***** error * sum of magnetic quantum number is not  ', &
!     & 'zero in cleb * m1=',f6.1,2x,'m2=',f6.1,2x,'m3=',f6.1,' *****')
!  601 format(' ***** error * sum of angular momentum is not integer ', &
!     & ' in cleb * j1=',f6.1,2x,'j2=',f6.1,2x,'j3=',f6.1,' *****')
!  602 format(' ***** error * triangular condition is biolated in cleb', &
!     & ' * j1=',f6.1,2x,'j2=',f6.1,2x,'j3=',f6.1,' *****')
!  603 format(' ***** error * (j,m) pairs invalid in cleb * j1=',f5.1, &
!     & 2x,'m1=',f6.1,2x,'j2=',f6.1,2x,'m2=',f6.1,2x,'j3=',f6.1,2x,    &
!     & 'm3=',f6.1,' *****')
!  605 format(' ***** error * abs(m.q.n.)  >  a.m. in cleb * j1=', &
!     & f6.1,', m1=',f6.1,', j1=',f6.1,', m2=',f6.1,', j3=',f6.1,     &
!     & ', m3 =',f6.1,' *****')
!  606 format(' ***** error * there are some negative a.m. in cleb * ', &
!     & 'j1=',f6.1,2x,'j2=',f6.1,2x,'j3=',f6.1,' *****')
 1000 return
end function clebsh
!=======================================================================



!=======================================================================
subroutine glweight(n,xmin,xmax,pmod,wmod)
      use consts,only:pi,epslim
!=======================================================================
! routine to get the points and weights of Gauss-Legendre quadrature
! coded by Kazuki Yoshida.
!   input:
!       n: number of points
!    xmin: minimum value of the integration range
!    xmax: maximum value of the integration range
! output:
!    pmod: node points
!    wmod: weight for each node point
!=======================================================================
      implicit real*8(a-h,o-z)
      integer,intent(in) :: n
      real*8,intent(in) :: xmin,xmax
      real*8,intent(out) :: pmod(n),wmod(n)
      real*8 p(n),w(n)
!=======================================================================
      imax=(n+1)/2
      do i=1,imax
! first candidate of ith root of Legendre polynomial p_{n}(x)=0.
       x=cos(pi*(i-0.250d0)/(n+0.50d0))
       do !loop for newton's method
!Legendre polynomial pl.
        pl1=1.0d0
        pl2=0.0d0
        do j=1,n
         pl3=pl2
         pl2=pl1
         pl1=((2.0d0*j-1.0d0)*x*pl2-(j-1.0d0)*pl3)/j
        end do
!derivative of Legendre polynomial dp.
        dp=n*(x*pl1-pl2)/(x**2.0d0-1.0d0)
        x1=x
        x=x1-pl1/dp
!convergence check for newton's method
        if(abs(x-x1)<epslim) exit
       end do
       p(i)=-x
       p(n+1-i)=x
       w(i)=2.0d0/((1.0d0-x**2.0d0)*dp**2.0d0)
       w(n+1-i)=w(i)
      end do
! changing the range; -1 to 1 -> xmin to xmax
      pmod(1:n)=(xmax-xmin)*p(1:n)/2.0d0+(xmax+xmin)/2.0d0
      wmod(1:n)=w(1:n)*(xmax-xmin)/2.d0

      return

end subroutine glweight
!=======================================================================



!=======================================================================
function polintm2d(x1a,x2a,ya,nx1,nx2,x1,x2,m1,m2)
!=======================================================================
! function to get m-point(m-1th order) polynomial interpolation
! for 2-dimension.
! coded by Kazuki Yoshida.
!  input:
!    nx1,nx2: number of points of x1_i and x2_j
!        x1a: array of points x1_i for i=1-n
!        x2a: array of points x2_j for j=1-n
!         ya: array of points f(x1_i,x2_j) for i,j=1-n
!      x1,x2: the point you want the value of f(x1,x2)
!      m1,m2: number of points that are used for polynomial interpolation
! output:
!  polintm2d: interpolated value f(x1,x2)
!=======================================================================
      implicit real*8(a-h,o-z)
      integer,intent(in) :: nx1,nx2,m1,m2
      real*8,intent(in) :: x1a(nx1),x2a(nx2),x1,x2,ya(nx1,nx2)
      real*8 y2m2(nx2)
!=======================================================================

      do ix2=1,nx2
       y2m2(ix2)=polintm(x1a,ya(:,ix2),nx1,x1,m1)
      end do
      polintm2d=polintm(x2a,y2m2,nx2,x2,m2)

      return

end function polintm2d
!=======================================================================



!=======================================================================
function polintm(xa,ya,n,x,m)
!=======================================================================
! function to get m-point(m-1th order) polynomial interpolation.
! coded by Kazuki Yoshida.
!   input:
!       n: number of points
!      xa: array of points x_i for i=1-n
!      ya: array of points f(x_i) for i=1-n
!       x: the point you want the value of f(x)
!       m: number of points that are used for polynomial interpolation
!          of f(x). x is in the middle of those points.
! output:
! polintm: interpolated value f(x)
!=======================================================================
      implicit real*8(a-h,o-z)
      integer,intent(in) :: n,m
      real*8,intent(in) :: x,xa(n),ya(n)
      real*8 xm(n),ym(n)
!=======================================================================
      if(m>n) then
        write(*,*) 'error in polintm: m > n'
        stop 1
      end if
! find j such that x is in between xa(j) and xa(j+1)
      call locate(xa,n,x,j)
! find k such that j is in the middle of k and k+m
      k=min(max(j-(m-1)/2,1),n+1-m)
      xm(1:m)=xa(k:k+m-1)
      ym(1:m)=ya(k:k+m-1)
      call polint(xm,ym,m,x,y)

      polintm=y

      return

end function polintm
!=======================================================================



!=======================================================================
subroutine polint(xa,ya,n,x,y)
!=======================================================================
! routine for polynomial interpolation with neville's algorithm,
! taken from numerical recipies,
! coded by Kazuki Yoshida.
!   input:
!       n: number of points
!      xa: array of points xa(i)=x_i for i=1-n
!      ya: array of points ya(i)=f(x_i) for i=1-n
!       x: the point where you want the value of f(x)
! output:
!       y: the value f(x)
!=======================================================================
      implicit real*8(a-h,o-z)
      integer,intent(in) :: n
      real*8,intent(in) :: xa(n),ya(n),x
      real*8,intent(out) :: y
      real*8 ca(n),da(n)
!=======================================================================
      nfound=1
      dif0=abs(x-xa(1))
      do i=1,n
       dif=abs(x-xa(i))
       if(dif<dif0) then
         nfound=i
         dif0=dif
       end if
       ca(i)=ya(i)
       da(i)=ya(i)
      end do
      y=ya(nfound)
      nfound=nfound-1
      do k=1,n-1
       do i=1,n-k
        xix=xa(i)-x
        xikx=xa(i+k)-x
        cddiff=ca(i+1)-da(i)
        xdiff=xix-xikx
        if(xdiff==0.0d0) then
          write(*,*) 'ERROR: xdiff=0 in polint'
          stop 1
        end if
        xdiff=cddiff/xdiff
        da(i)=xikx*xdiff
        ca(i)=xix*xdiff
       end do
       if(2*nfound<n-k) then
         dy=ca(nfound+1)
        else
         dy=da(nfound)
         nfound=nfound-1
       end if
       y=y+dy
      end do

      return

end subroutine polint
!=======================================================================



!=======================================================================
subroutine locate(xx,n,x,j)
!=======================================================================
! routine to find index j such that x is between xx(j) and xx(j+1),
! taken from numerical recipes,
! coded by Kazuki Yoshida.
!   input:
!      xx: array of points xx(i)=x_i for i=1-n
!       n: number of points
!       x: the point between xx(j) and xx(j+1)
! output:
!       j: the index of the array xx such that xx(j) < x < xx(j+1)
!=======================================================================
      implicit real*8(a-h,o-z)
      real*8 xx(n)
!=======================================================================
      ji=0
      jf=n+1
      do
       if(jf-ji<=1) then
         exit
        else
         jcal=(jf+ji)/2
         if((xx(n)>=xx(1)).eqv.(x>=xx(jcal))) then
           ji=jcal
          else
           jf=jcal
         end if
       end if
      end do
      if(x==xx(1)) then
        j=1
       else
        j=ji
      end if

      return

end subroutine locate
!=======================================================================



!=======================================================================
!**
!**        ***** suphod2d *****
!**
!    this program is o.k. only when jisu=3
!                                    and dx1,dx2=const  -- care --

function suphod2d(x1in,x1,dx1,x2in,x2,dx2,y,n1max,n2max,n1dim,n2dim)
      use kibmod
      implicit real*8(a-h,o-z)
      real*8 x1(n1dim),x2(n2dim),y(n1dim,n2dim)
      integer,parameter :: jisu=3
!** ===========================================
      if(n1max<=jisu .or. n2max<=jisu) then
        write(kibout,400) n1max,n2max,jisu
  400   format(/1x,' n1max,n2max,jisu=',3i5,10x,'stop in suphod2d')
        stop 1
      end if

      nnr1=int((x1in-x1(1))/dx1+1.01)
      if(nnr1<=0 .or. nnr1>n1max) then
        write(kibout,500) nnr1
  500   format(/1x,' nnr1=',i5,10x,'stop in suphod2d --- x1in < x1(1) ' &
     &            ,' or  x1in > x1(n1max) ---')
        stop 1
      end if

      nnr2=int((x2in-x2(1))/dx2+1.01)
      if(nnr2<=0 .or. nnr2>n2max) then
        write(*,*) x2in,x2(1),dx2,x2(n2max)
        write(kibout,501) nnr2
  501   format(/1x,' nnr2=',i5,10x,'stop in suphod2d --- x2in < x2(1) ' &
     &            ,' or  x2in > x2(n2max) ---')
        stop 1
      end if

      j1min=nnr1-1
      j1max=nnr1+2
      if(j1min<=1) j1min=1
      if(j1max>=n1max) j1min=n1max-3

      x10=x1(j1min)
      x11=x1(j1min+1)
      x12=x1(j1min+2)
      x13=x1(j1min+3)

      j2min=nnr2-1
      j2max=nnr2+2
      if(j2min<=1) j2min=1
      if(j2max>=n2max) j2min=n2max-3

      x20=x2(j2min)
      x21=x2(j2min+1)
      x22=x2(j2min+2)
      x23=x2(j2min+3)

      x1in0=x1in-x10
      x1in1=x1in-x11
      x1in2=x1in-x12
      x1in3=x1in-x13

      x2in0=x2in-x20
      x2in1=x2in-x21
      x2in2=x2in-x22
      x2in3=x2in-x23

      x101=x10-x11
      x102=x10-x12
      x103=x10-x13

      x110=x11-x10
      x112=x11-x12
      x113=x11-x13

      x120=x12-x10
      x121=x12-x11
      x123=x12-x13

      x130=x13-x10
      x131=x13-x11
      x132=x13-x12

      x201=x20-x21
      x202=x20-x22
      x203=x20-x23

      x210=x21-x20
      x212=x21-x22
      x213=x21-x23

      x220=x22-x20
      x221=x22-x21
      x223=x22-x23

      x230=x23-x20
      x231=x23-x21
      x232=x23-x22

      y00=y(j1min,j2min)
      y01=y(j1min,j2min+1)
      y02=y(j1min,j2min+2)
      y03=y(j1min,j2min+3)

      x2x0= x2in1/x201 *x2in2/x202 *x2in3/x203
      x2x1= x2in0/x210 *x2in2/x212 *x2in3/x213
      x2x2= x2in0/x220 *x2in1/x221 *x2in3/x223
      x2x3= x2in0/x230 *x2in1/x231 *x2in2/x232

      y0=x2x0*y00 +x2x1*y01 +x2x2*y02 +x2x3*y03

      y10=y(j1min+1,j2min)
      y11=y(j1min+1,j2min+1)
      y12=y(j1min+1,j2min+2)
      y13=y(j1min+1,j2min+3)

      y1=x2x0*y10 +x2x1*y11 +x2x2*y12 +x2x3*y13

      y20=y(j1min+2,j2min)
      y21=y(j1min+2,j2min+1)
      y22=y(j1min+2,j2min+2)
      y23=y(j1min+2,j2min+3)

      y2=x2x0*y20 +x2x1*y21 +x2x2*y22 +x2x3*y23

      y30=y(j1min+3,j2min)
      y31=y(j1min+3,j2min+1)
      y32=y(j1min+3,j2min+2)
      y33=y(j1min+3,j2min+3)

      y3=x2x0*y30 +x2x1*y31 +x2x2*y32 +x2x3*y33

      x1x0= x1in1/x101 *x1in2/x102 *x1in3/x103
      x1x1= x1in0/x110 *x1in2/x112 *x1in3/x113
      x1x2= x1in0/x120 *x1in1/x121 *x1in3/x123
      x1x3= x1in0/x130 *x1in1/x131 *x1in2/x132


      suphod2d=x1x0*y0 +x1x1*y1 +x1x2*y2 +x1x3*y3

      return

end function suphod2d
!=======================================================================



!=======================================================================
subroutine numerovl(ll,gx,dx,nxmax,y)
      implicit real*8 (a-h,o-z)
      integer,intent(in) :: ll,nxmax
      real*8, intent(in) :: dx
      complex*16,intent(in) :: gx(nxmax)
      complex*16,intent(out) :: y(nxmax)
!=======================================================================
!      routine to solve the differential equation
!      d^2y/dx^2 = g(x)y(x)
!      with numerov's method. coded by K.Yoshida.
!
!      input:
!          ll: angular momentum of the partial wave
!          gx: g(x) shown above
!          dx: bin size
!       nxmax: size of array
!     output:
!           y: y(x)
!=======================================================================

      y(:)=(0.d0,0.d0)
      y(2)=dcmplx(dx**(ll+1), 0.d0)

      dx2=dx**2
      dx256=dx2*5.0d0/6.0d0
      dx212=dx2/12.0d0

!special condition for l=1 at x=0 due to 1/r**2 behavior of l barrier.
      if(ll==1) then
        y(3)=(2.d0+dx256*gx(2))*y(2) &
     &      +dx212*2.d0*y(2)/(dx**(ll+1))
        y(3)=y(3)/(1.d0-dx212*gx(3))
       else
        y(3)=(2.d0+dx256*gx(2))*y(2) &
     &      -(1.d0-dx212*gx(1))*y(1)
        y(3)=y(3)/(1.d0-dx212*gx(3))
      end if

      do ix=4,nxmax
       y(ix)=(2.d0+dx256*gx(ix-1))*y(ix-1) &
     &      -(1.d0-dx212*gx(ix-2))*y(ix-2)
       y(ix)=y(ix)/(1.d0-dx212*gx(ix))
      end do

      return

end subroutine numerovl
!=======================================================================



!=======================================================================
subroutine potint(nrmax,nepfmax,iepfix,iedg,epfmin,epfix,tin,detbl,beta &
     &           ,uopttbl,ulstbl,fnltbl,uopt,uoptls,fnl)
      use consts
      implicit real*8 (a-h,o-z)
      integer,intent(in) :: nrmax,nepfmax,iepfix,iedg
      real*8,intent(in) :: epfmin,epfix,tin,detbl,beta
      complex*16,intent(in) ::  uopttbl(nepfmax,nrmax) &
     &                         ,ulstbl(nepfmax,nrmax),fnltbl(nepfmax,nrmax)
      complex*16,intent(out) :: uopt(nrmax),uoptls(nrmax),fnl(nrmax)
!=======================================================================
      if(iepfix==1 .and. abs(epfix-tin)>epslim .and. iedg==0) then
        write(*,'(a,f0.4,a)') 'ERROR: potential at energy', tin, ' not given'
        stop 1
      end if
      if(iepfix/=1) then
        iecal=int((tin-epfmin)/detbl+1.01)
        if(iecal<1 .or. iecal>=nepfmax) then
          if(iedg==0 .and. iecal/=nepfmax) then
            write(*,'(a,f0.4,a)') 'ERROR: potential at energy', tin, ' not given'
            stop 1
          end if
          if(iecal<1) iec=1
          if(iecal>=nepfmax) iec=nepfmax
          do ir=1,nrmax
           uopt(ir)=uopttbl(iec,ir)
           uoptls(ir)=ulstbl(iec,ir)
           if(beta<-epslim) fnl(ir)=fnltbl(iec,ir)
          end do
         else
          ebs=(iecal-1)*detbl+epfmin
          do ir=1,nrmax
           uopt(ir)=uopttbl(iecal,ir)                         &
     &               +(uopttbl(iecal+1,ir)-uopttbl(iecal,ir)) &
     &                *(tin-ebs)/detbl
           uoptls(ir)=ulstbl(iecal,ir)                        &
     &                 +(ulstbl(iecal+1,ir)-ulstbl(iecal,ir)) &
     &                 *(tin-ebs)/detbl
           if(beta<-epslim)                                   &
     &       fnl(ir)=fnltbl(iecal,ir)                         &
     &                 +(fnltbl(iecal+1,ir)-fnltbl(iecal,ir)) &
     &                  *(tin-ebs)/detbl
          end do
        end if
      end if

      return

end subroutine potint
!=======================================================================



!=======================================================================
subroutine coulfg(l,eta,rho1,drho,f1,g1,fp1,gp1)
      implicit real*8(a-h,o-z)

      real*8 f(500),g(500),fp(500),gp(500)

      drhofg=dmin1(0.2d0,drho)
      rho2=rho1+drhofg

      ll=l+1 +1
      e=eta
      h=drhofg
      r=rho1
      rp=rho2
      te=e+e
      tf=e**2
      if(ll-50) 20,35,35
   20 elp=50.
      j=50
      go to 45
   35 elp=ll
      j=ll
   45 a=datan(e/elp)
      b=dsqrt(tf+elp**2)
      y=a*(elp-0.5d0)+e*(dlog(b)-1.0d0)-dsin(a)/(12.0d0*b)          &
     &  +dsin(3.0d0*a)/(360.0d0*b**3)-dsin(5.0d0*a)/(1260.0d0*b**5) &
     &  +dsin(7.0d0*a)/(1680.0d0*b**7)-dsin(9.0d0*a)/(1188.0d0*b**9)
      k=j-1
      if(j-ll)65,65,70
   65 s1=y
   70 do 100 i=1,k
      elp=elp-1.
      j=j-1
      y=y-datan(e/elp)
  100 continue
      s1=y
      del1=r-te
      rmax=dmax1(10.0d0,(tf+3.0d0+4.0d0*e)*5.0d0/12.0d0)
      del=r-rmax
      if(e-5.)280,130,130
  130 if(dabs(del1)-dabs(del))140,140,280
  140 del=del1
      if(del)147,145,147
  145 i=2
      go to 150
  147 i=1
  150 x=te
      t1=tf
      t2=t1**2
!**   t3=e** .666666667
!**   t9=e** .166666667
      t9=e**(1.0d0/6.0d0)
      t3=t9**4
      t4=t3**2
      t5=t4**2
      t6=t3*t5
      t7=t4*t6
      t8=t3*t7
      y=1.22340402d0*t9*(1.0d0+0.495957017d-1/t4- &
     &  0.888888889d-2/t1+0.245519918d-2          &
     &  /t6-0.910895806d-3/t2+0.253468412d-3/t8)
      z=-0.707881773d0/t9*(1.0d0-0.172826039d0/t3+0.317460317d-3/t1- &
     &  0.358121485d-2/t5+0.311782468d-3/t2-0.907396643d-3/t7)
      go to 665
  280 if(e)285,290,285
  285 if(del)310,290,290
  290 x=r
      i=2
      go to 320
  310 x=rmax
      i=1
  320 t1=tf
      t2=x+x
      t3=x-e*dlog(t2)+s1
      t4=e/t2
      ss=1.0d0
      ts=0.0
      sl=0.0
      tl=1.0d0-e/x
      sss=1.0d0
      sts=0.0
      ssl=0.0
      stl=tl
      en=0.0
      do 620 k=1,25
      t5=en+1.
      t6=t5+en
      t7=en*t5
      t8=t6*t4/t5
      t9=(t1-t7)/(t2*t5)
      t5=t8*ss-t9*ts
      ts=t8*ts+t9*ss
      ss=t5
      if(dabs(ss/sss)-1.0d-10) 630,630,540
  540 t5=t8*sl-t9*tl-ss/x
      tl=t8*tl+t9*sl-ts/x
      sl=t5
      sss=sss+ss
      sts=sts+ts
      ssl=ssl+sl
      stl=stl+tl
      en=en+1.
  620 continue
  630 t8=dsin(t3)
      t9=dcos(t3)
      y=sss*t9-sts*t8
      z=ssl*t9-stl*t8
  665 go to (670,810),i
  670 m=1
  671 n=nint(dabs(del/h))
      if(n)675,675,680
  675 dx=del
      go to 700
  680 en=n
      dx=del/en
  700 t1=0.5d0*dx
      t2=0.25d0*t1
      t3=te
      do 805 i=1,n
      t4=dx*(t3/x-1.)*y
      x=x+t1
      t5=dx*(t3/x-1.)*(y+t1*z+t2*t4)
      x=x+t1
      t6=dx*(t3/x-1.)*(y+dx*z+t1*t5)
      y=y+dx*(z+(t4+t5+t5)/6.0d0)
      z=z+(t4+4.0d0*t5+t6)/6.0d0
  805 continue
      go to (810,828),m
  810 g(1)=y
      m=2
      del=rp-r
      w=z
      go to 671
  828 gp(1)=y
      t1=tf
      t8=dsqrt(1.+t1)
      g(2)=((1./r+e)*g(1)-w)/t8
      gp(2)=((1./rp+e)*y-z)/t8
      t2=1.0d0
      t3=2.0d0
      do 910 i=3,ll
      t4=t2+t3
      t5=t2*t3
      t6=t3*dsqrt(t2**2+t1)
      t7=t2*dsqrt(t3**2+t1)
      g (i)=(t4*(e+t5/r )*g (i-1)-t6*g (i-2))/t7
      gp(i)=(t4*(e+t5/rp)*gp(i-1)-t6*gp(i-2))/t7
      t2=t2+1.0d0
      t3=t3+1.0d0
  910 continue
!**   i=l+11
      i=ll+11
      n=nint(r+r+11.0d0)
      if(i-n)960,960,950
  950 n=i
  960 y=1.0d-20
      yp=y
      x=y
      xp=x
      z=0.0
      zp=z
      t2=n
 1000 t3=t2+1.0d0
      t4=t2+t3
      t5=t2*t3
      t6=t2*dsqrt(t3**2+t1)
      t7=t3*dsqrt(t2**2+t1)
      y =(t4*(e+t5/r )*y -t6*z )/t7
      yp=(t4*(e+t5/rp)*yp-t6*zp)/t7
      if(n-ll)1060,1060,1080
 1060 f(n)=y
      fp(n)=yp
      go to 1120
 1080 if(1.0d0-dabs(y))1090,1090,1120
 1090 y=y*1.0d-20
      yp=yp*1.0d-20
      x=x*1.0d-20
      xp=xp*1.0d-20
 1120 n=n-1
      z=x
      zp=xp
      x=y
      xp=yp
      t2=t2-1.0d0
      if(n)1150,1150,1000
 1150 y=f(1)*g(2)-f(2)*g(1)
      yp=fp(1)*gp(2)-fp(2)*gp(1)
      z=1.0d0/(y*t8)
      zp=1.0d0/(yp*t8)
      do 1180 i=1,ll
      fp(i)=fp(i)*zp
 1180 f(i)=f(i)*z
      f1=f(l+1)

!----f2 and g2 are not used.
!     f2=fp(l+1)
      g1=g(l+1)
!     g2=gp(l+1)
      al1=l+1
      a=al1**2/rho1+eta
      b=al1**2+eta**2
      b=dsqrt(b)
      fp1=(a*f1-b*f(l+2))/al1
      gp1=(a*g1-b*g(l+2))/al1
      return

end subroutine coulfg
!=======================================================================



!=======================================================================
function dfunc(a,b,c,cost)
      use consts,only:pi,epslim
!=======================================================================

!      reduced wiger-d matrix calculation code
!      written by K. Yoshida.
!
!      cost=cos(tht),
!      dfunc(a,b,c,cost) = d^{a}_{bc}(tht)
!
!      calculation has been confirmed individually for
!      a=1/2, b,c=-1/2,1/2
!      a=1  , b,c=-1,0,1
!
!      orthogonal relations
!      sum_{k} d^{l}_{mk}(tht) d^{l}_{m'k}(tht) = \delta_{mm'}
!      sum_{m} d^{l}_{mk}(tht) d^{l}_{mk'}(tht) = \delta_{kk'}
!      are confirmed for integer and half-integer l,m,k < 75.
!
!      symmetric property
!      d^{a}_{bc}(tht)=(-)^{m-k}d^{a}_{cb}(tht)
!      is confirmed up to integer and half-integer l,m,k < 75.

!=======================================================================
      implicit real*8 (a-h,o-z)
      integer,parameter :: lmax=100
      real*8, intent(in) :: a,b,c,cost
      real*8 dlmk0(-2:2,-2:2)
      real*8 dlmk(-2*lmax-4:2*lmax+4,0:1)
      real*8 redcgm1(0:1),redcgk(-1:1)
!=======================================================================

      dfunc=0.d0

      if(a>dble(lmax)) then
        write(*,*) 'ERROR: a > lmax in function dfunc'
        stop 1
      end if

      if(abs(cost)>1.0d0+epslim) then
        write(*,*) 'ERROR: abs(cost)>1.d0'
        stop 1
      end if
      dlmk0(:,:)=0.d0
      dlmk(:,:)=0.d0

      na2=nint(2.d0*a)
      nb2=nint(2.d0*b)
      nc2=nint(2.d0*c)

! input values check
      if(abs(na2-2.d0*a)>epslim) then
        write(*,*) 'ERROR: a value is not integer nor half-integer'
        stop 1
      end if

      if(abs(nb2-2.d0*b)>epslim) then
        write(*,*) 'ERROR: b value is not integer nor half-integer'
        stop 1
      end if

      if(abs(nc2-2.d0*c)>epslim) then
        write(*,*) 'ERROR: c value is not integer nor half-integer'
        stop 1
      end if

      if(na2<1) then
        if(na2==0) then
          dfunc=1.d0
          return
         else
          write(*,*) na2
          write(*,*) 'ERROR: in the integer a'
          stop 1
        end if
      end if

      if(abs(nb2)>na2) then
        write(*,*) 'ERROR: b is larger than a'
        stop 1
      end if

      if(abs(nc2)>na2) then
        write(*,*) 'ERROR: c is larger than a'
        stop 1
      end if

      if(mod(na2+nb2,2)==1) then
        write(*,*)'ERROR: b should be (half-)integer when a is (half-)integer'
        stop 1
      end if

      if(mod(na2+nc2,2)==1) then
        write(*,*)'ERROR: c should be (half-)integer when a is (half-)integer'
        stop 1
      end if
! input values check end

      csh=sqrt((1.d0+cost)/2.d0) !cos(tht/2)
      snh=sqrt((1.d0-cost)/2.d0) !sin(tht/2)

      dlmk0(-1,-1)= csh
      dlmk0(-1, 1)= snh
      dlmk0( 1,-1)=-snh
      dlmk0( 1, 1)= csh

      if(na2==1) then
        dfunc=dlmk0(nb2,nc2)
        return
      end if

! reduction of the calculation using the symmetry,
! and the determination of the factor.
!---------------------------------------------- begin
      cosdum=cost
      if(nb2>=0) then
        if(nc2>=0) then
          fac=1.d0
         else
          nc2=-nc2
          cosdum=-cost
          fac=(-1)**((na2+nb2)/2)
        end if
       else
        if(nc2>0) then
          nb2=-nb2
          cosdum=-cost
          fac=(-1)**((na2-nc2)/2)
         else
          nb2=-nb2
          nc2=-nc2
          fac=(-1)**((nb2-nc2)/2)
        end if
      end if

      if(nc2>=nb2) then
        nbcal2=nb2
        nccal2=nc2
        fac2=fac
       else
        nbcal2=nc2
        nccal2=nb2
        fac2=fac*(-1)**((nb2-nc2)/2)
      end if
!---------------------------------------------- end

      csh=sqrt((1.d0+cosdum)/2.d0)
      snh=sqrt((1.d0-cosdum)/2.d0)

      dlmk0( 0,-2)=-sqrt(2.d0)*csh*snh
      dlmk0( 0, 0)=-snh**2.d0+csh**2.d0
      dlmk0( 0, 2)= sqrt(2.d0)*csh*snh
      dlmk0( 2,-2)= snh**2.d0
      dlmk0( 2, 0)=-sqrt(2.d0)*csh*snh
      dlmk0( 2, 2)= csh**2.d0

      if(mod(na2,2)==1) then  ! initial value for the loop
        l0=3
       else
        l0=2
      end if

      if(l0==2) then
        dlmk(0,0)=1.d0
       else
        dlmk(-1,0)=-snh
        dlmk( 1,0)= csh
      end if

      kc=0
      kd=1
      mmin=l0-2

      do l2=l0,na2,2 ! find dlmk recursively
       mdum=mmin
       mmin=min(l2,nbcal2)
       m1=(mmin-mdum)/2
       fl=l2/2.d0-1.d0
       fm=mmin/2.d0
       fdeno=sqrt((2.d0*fl+1.d0)*(2.d0*fl+2.d0))
       redcgm1(0)=sqrt((fl-fm+1.d0)*(fl+fm+1.d0)*2.d0)/fdeno
       redcgm1(1)=sqrt((fl+fm)*(fl+fm+1.d0))/fdeno

       kmin=max(-l2,nccal2-na2+l2)
       kmax=min(l2,nccal2+na2-l2)

       do ki=kmin,kmax,2
        fk=ki/2.d0
        redcgk(-1)=sqrt((fl-fk)*(fl-fk+1.d0))/fdeno
        redcgk( 0)=sqrt((fl-fk+1.d0)*(fl+fk+1.d0)*2.d0)/fdeno
        redcgk( 1)=sqrt((fl+fk)*(fl+fk+1.d0))/fdeno

        fsum=0.d0

        do k=-1,1
         fsum=fsum+redcgk(k)*dlmk0(m1*2,k*2)*dlmk(ki-k*2,kc)
         dlmk(ki,kd)=fsum/redcgm1(m1)
        end do

       end do

       dlmk(kmin-2,kd)=0.d0
       dlmk(kmin-4,kd)=0.d0
       dlmk(kmax+2,kd)=0.d0
       dlmk(kmax+4,kd)=0.d0

       ke=kd
       kd=kc
       kc=ke

      end do

      dfunc=dlmk(nccal2,kc)*fac2

      return

end function dfunc
!=======================================================================



!=======================================================================
subroutine nnkincal(q0,q,thq0,thq,phq0,phq,ielmedg,ionsh ,excd &
     &             ,nq0cal,nqcal,nphq0cal,nphqcal              &
     &             ,q0cal,qcal,thq0cal,thqcal,phq0cal,phqcal)

      use nntbl,only:q0mn,q0mx,dq0,qmn,qmx,dq &
     &             ,thq0mn,thq0mx,thqmn,thqmx &
     &             ,phq0mn,phq0mx,dphq0,phqmn,phqmx,dphq
      use kibmod,only:kibout
!------------------------------------------------------------------------------------
      implicit none
      real*8,intent(in) :: q0,q,thq0,thq,phq0,phq
      integer, intent(in) :: ielmedg,ionsh
      integer, intent(out) :: nq0cal,nqcal,nphq0cal,nphqcal
      real*8,intent(out) :: q0cal,qcal,thq0cal,thqcal,phq0cal,phqcal
      integer ilp
      character*1 excd(6)
!------------------------------------------------------------------------------------

      q0cal=q0
      qcal=q
      thq0cal=thq0
      thqcal=thq
      phq0cal=phq0
      phqcal=phq

      excd(:)=' '

      if(q<qmn .or. q>qmx) then
        excd(1)='%'
        if(ionsh==0 .or. ionsh==1) excd(1)='*'
      end if
      if(q0<q0mn .or. q0>q0mx) then
        excd(2)='%'
        if(ionsh==0 .or. ionsh==2) excd(2)='*'
      end if
      if(thq<thqmn .or. thq>thqmx)     excd(3)='*'
      if(thq0<thq0mn .or. thq0>thq0mx) excd(4)='*'
      if(phq<phqmn .or. phq>phqmx)     excd(5)='*'
      if(phq0<phq0mn .or. phq0>phq0mx) excd(6)='*'

      if(q0<q0mn) q0cal=q0mn
      if(q0>q0mx) q0cal=q0mx
      if(q<qmn) qcal=qmn
      if(q>qmx) qcal=qmx

      if(thq0<thq0mn) thq0cal=thq0mn
      if(thq0>thq0mx) thq0cal=thq0mx
      if(thq<thqmn) thqcal=thqmn
      if(thq>thqmx) thqcal=thqmx

      if(phq0<phq0mn) phq0cal=phq0mn
      if(phq0>phq0mx) phq0cal=phq0mx
      if(phq<phqmn) phqcal=phqmn
      if(phq>phqmx) phqcal=phqmx

      nq0cal  =nint((q0cal-q0mn)/dq0)+1
      nqcal   =nint((qcal-qmn)/dq)+1
      nphq0cal=nint((phq0cal-phq0mn)/dphq0)+1
      nphqcal =nint((phqcal-phqmn)/dphq)+1

      if(ielmedg==0) then
        do ilp=1,6
         if(excd(ilp)=='*') then
           write(kibout,601) q,qmn,qmx,q0,q0mn,q0mx             &
     &                      ,thq,thqmn,thqmx,thq0,thq0mn,thq0mx &
     &                      ,phq,phqmn,phqmx,phq0,phq0mn,phq0mx
  601      format('ERROR: kinematical variable(s) not prepared in the table' &
     &           /'q    qmn    qmx   :',3f10.4                               &
     &           /'q0   q0mn   q0mx  :',3f10.4                               &
     &           /'thq  thqmn  thqmx :',3f10.4                               &
     &           /'phq  phqmn  phqmx :',3f10.4                               &
     &           /'thq0 thq0mn thq0mx:',3f10.4                               &
     &           /'phq0 phq0mn phq0mx:',3f10.4)
           stop 1
         end if
        end do
      end if

end subroutine nnkincal
!=======================================================================



!=======================================================================
subroutine makenntbl()
!     This routine will be called only once but returned values should be kept
!     during entire run. All arguments of this routine is stored in module nntbl.
      use dims,only:nq0dim,nthq0dim,nthqdim,nphq0dim,nphqdim
      use consts,only:pi
      use nntbl

!------------------------------------------------------------------------------------
      implicit none
      integer i
!------------------------------------------------------------------------------------

      nq0mn=1
      nq0mx=nint((q0mx-q0mn)/dq0)+1
      nqmn=1
      nqmx=nint((qmx-qmn)/dq)+1

      thq0mn=th0mn*pi/180.d0
      thq0mx=th0mx*pi/180.d0
      dthq0=dth*pi/180.d0
      phq0mn=ph0mn*pi/180.d0
      phq0mx=ph0mx*pi/180.d0
      dphq0=dph*pi/180.d0
      nth0mn=1
      nph0mn=1
      nth0mx=nint((thq0mx-thq0mn)/dthq0)+1
      nph0mx=nint((phq0mx-phq0mn)/dphq0)+1

      thqmn=thmn*pi/180.d0
      thqmx=thmx*pi/180.d0
      dthq=dth*pi/180.d0
      phqmn=phmn*pi/180.d0
      phqmx=phmx*pi/180.d0
      dphq=dph*pi/180.d0
      nthmn=1
      nphmn=1
      nthmx=nint((thqmx-thqmn)/dthq)+1
      nphmx=nint((phqmx-phqmn)/dphq)+1

      if(nq0mx>nq0dim) then
        write(*,*) 'ERROR: nq0mx > nq0dim'
        write(*,'(a,i0,a,i0)') ' nq0mx=',nq0mx,', nq0dim=',nq0dim
        stop 1
      end if
      if(nth0mx>nthq0dim) then
        write(*,*) 'ERROR: nth0mx > nthq0dim'
        write(*,'(a,i0,a,i0)') ' nth0mx=',nth0mx,', nthq0dim=',nthq0dim
        stop 1
      end if
      if(nthmx>nthqdim) then
        write(*,*) 'ERROR: nthmx > nthqdim'
        write(*,'(a,i0,a,i0)') ' nthmx=',nthmx,', nthqdim=',nthqdim
        stop 1
      end if
      if(nph0mx>nphq0dim) then
        write(*,*) 'ERROR: nph0mx > nphq0dim'
        write(*,'(a,i0,a,i0)') ' nph0mx=',nphmx,', nphq0dim=',nphq0dim
        stop 1
      end if
      if(nphmx>nphqdim) then
        write(*,*) 'ERROR: nphmx > nphqdim'
        write(*,'(a,i0,a,i0)') ' nphmx=',nphmx,', nphqdim=',nphqdim
        stop 1
      end if

      do i=nthmn,nthmx
       thqtbl(i)=(i-nthmn)*dthq+thqmn
      end do
      do i=nth0mn,nth0mx
       thq0tbl(i)=(i-nth0mn)*dthq+thq0mn
      end do
      do i=nphmn,nphmx
       phqtbl(i)=(i-nphmn)*dphq+phqmn
      end do
      do i=nph0mn,nph0mx
       phq0tbl(i)=(i-nph0mn)*dphq+phq0mn
      end do

end subroutine makenntbl
!=======================================================================



!=======================================================================
subroutine getelmnntbl(nphqcal,nq0cal,nphq0cal,thqcal,thq0cal,elmnntbl)
! Make a table of t_{NN} for each spin configuration(ispin).
! Interpolation for theta(thq0cal) and theta'(thqcal) is done in this routine.

      use array,only:gma
      use nntbl
!------------------------------------------------------------------------------------
      implicit none
      integer,intent(in) :: nphqcal,nq0cal,nphq0cal
      real*8,intent(in) :: thqcal,thq0cal
      complex*16,intent(out) :: elmnntbl(16)
      integer iupsp,iup1fz,iup2fz,iup3fz,ispin
      real*8 polintm2d
!------------------------------------------------------------------------------------

      do iupsp=1,2
       do iup1fz=1,2
        do iup2fz=1,2
         do iup3fz=1,2

          ispin=(iup1fz-1)*8+(iupsp-1)*4+(iup2fz-1)*2+(iup3fz-1)*1+1

          elmnntbl(ispin)=                                              &
     &     dcmplx(polintm2d(thqtbl(nthmn:nthmx),thq0tbl(nth0mn:nth0mx)  &
     &                     ,dble(gma(nthmn:nthmx,nphqcal                &
     &                          ,nq0cal,nth0mn:nth0mx,nphq0cal,ispin))  &
     &                     ,nthmx,nth0mx,thqcal,thq0cal,3,3)            &
     &           ,polintm2d(thqtbl(nthmn:nthmx),thq0tbl(nth0mn:nth0mx)  &
     &                     ,dimag(gma(nthmn:nthmx,nphqcal               &
     &                           ,nq0cal,nth0mn:nth0mx,nphq0cal,ispin)) &
     &                     ,nthmx,nth0mx,thqcal,thq0cal,3,3))

         end do
        end do
       end do
      end do

end subroutine getelmnntbl
!=======================================================================


!=======================================================================
subroutine header(ivar,ifrm,imir,ical,ixunt,kunt)
      use kibmod,only:kibtbl,kiblg,kibpx,kibtr,kibtl
      use dims
      implicit real*8 (a-h,o-z)
      character kpb1*2,kpb2*7,tqdx*20,mdu*13,mdu2*15,mdu3*15
!=======================================================================

      if(kunt==0) then
        kpb1='kb'
        kpb2='[1/fm] '
       else if(kunt==1) then
        kpb1='pb'
        kpb2='[MeV/c]'
       else if(kunt==2) then
        kpb1='pb'
        kpb2='[GeV/c]'
      end if

      if(ivar==1 .or. (ivar>9 .and. ivar<20) .or. ivar>39) then
        if(ixunt==1) then
          tqdx='tdx[ub/(MeVsr2)]    '
         else
          tqdx='tdx[mb/(MeVsr2)]    '
        end if
       else if(ivar==2 .or. (ivar>19 .and. ivar<30)) then
        if(ixunt==1) then
          if(kunt==0) then
            tqdx='tdx[ub/(fm-3sr2)]   '
           else if(kunt==1) then
            tqdx='tdx[ub/(MeV3c-3sr2)]'
           else if(kunt==2) then
            tqdx='tdx[ub/(GeV3c-3sr2)]'
          end if
         else
          if(kunt==0) then
            tqdx='tdx[mb/(fm-3sr2)]   '
           else if(kunt==1) then
            tqdx='tdx[mb/(MeV3c-3sr2)]'
           else if(kunt==2) then
            tqdx='tdx[mb/(GeV3c-3sr2)]'
          end if
        end if
       else if(ivar==3 .or. (ivar>29 .and. ivar<40)) then
        if(ixunt==1) then
          tqdx='qdx[ub/(MeV2srrad)] '
         else
          tqdx='qdx[mb/(MeV2srrad)] '
        end if
       else if(ivar==9) then
        if(ixunt==1) then
          if(kunt==0) then
            mdu='[ub/fm-1]    '
            mdu2='[ub/fm-2]      '
            mdu3='[ub/fm-3]      '
           else if(kunt==1) then
            mdu='[ub/(MeVc-1)]'
            mdu2='[ub/(MeVc-1)^2]'
            mdu3='[ub/(MeVc-1)^3]'
           else if(kunt==2) then
            mdu='[ub/(GeVc-1)]'
            mdu2='[ub/(GeVc-1)^2]'
            mdu3='[ub/(GeVc-1)^3]'
          end if
         else
          if(kunt==0) then
            mdu='[mb/fm-1]    '
            mdu2='[mb/fm-2]      '
            mdu3='[mb/fm-3]      '
           else if(kunt==1) then
            mdu='[mb/(MeVc-1)]'
            mdu2='[mb/(MeVc-1)^2]'
            mdu3='[mb/(MeVc-1)^3]'
           else if(kunt==2) then
            mdu='[mb/(GeVc-1)]'
            mdu2='[mb/(GeVc-1)^2]'
            mdu3='[mb/(GeVc-1)^3]'
          end if
        end if
      end if

      if(ivar/=9) then
        if(ifrm==0) then
          if(imir==0) then
            if(ical==0) then
              write(kibtbl,610) kpb1,kpb2
             else
              write(kibtbl,611) kpb1,kpb2,tqdx
            end if
           else
            if(ical==0) then
              write(kibtbl,612) kpb1,kpb2
             else
              write(kibtbl,613) kpb1,kpb2,tqdx
            end if
          end if
         else if(ifrm==1) then
          if(imir==0) then
            if(ical==0) then
              write(kibtbl,614) kpb1,kpb2
             else
              write(kibtbl,615) kpb1,kpb2,tqdx
            end if
           else
            if(ical==0) then
              write(kibtbl,616) kpb1,kpb2
             else
              write(kibtbl,617) kpb1,kpb2,tqdx
            end if
          end if
         else if(ifrm==2) then
          if(imir==0) then
            if(ical==0) then
              write(kibtbl,618) kpb1,kpb2
             else
              write(kibtbl,619) kpb1,kpb2,tqdx
            end if
           else
            if(ical==0) then
              write(kibtbl,620) kpb1,kpb2
             else
              write(kibtbl,621) kpb1,kpb2,tqdx
            end if
          end if
        end if
       else
        if(kiblg>0) then
          if(imir==0) then
            write(kiblg,701) kpb1,kpb2,mdu
           else
            write(kiblg,702) kpb1,kpb2,mdu
          end if
        end if
        if(kibpx>0) write(kibpx,703) kpb1,kpb2,mdu
        if(kibtr>0) write(kibtr,704) kpb1,kpb2,mdu,kpb1,mdu2
        if(kibtl>0) write(kibtl,705) kpb1,kpb2,mdu,kpb1,mdu3
      end if

  610 format(2x,'t1l[MeV]',3x,'th1l[deg]',3x,'ph1l[deg]',1x,'t2l[MeV]' &
     &      ,3x,'th2l[deg]',3x,'ph2l[deg]',1x,a2,'l',a7,1x,'thbl[deg]' &
     &      ,3x,'phbl[deg]',1x,'pr[MeV/c]',1x,'isol')
  611 format(2x,'t1l[MeV]',3x,'th1l[deg]',3x,'ph1l[deg]',1x,'t2l[MeV]' &
     &      ,3x,'th2l[deg]',3x,'ph2l[deg]',1x,a2,'l',a7,1x,'thbl[deg]' &
     &      ,3x,'phbl[deg]',1x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  612 format(2x,'t1l[MeV]',3x,'th1l:m[deg]',1x,'ph1l[deg]',1x,'t2l[MeV]'   &
     &      ,3x,'th2l:m[MeV]',1x,'ph2l[deg]',1x,a2,'l',a7,1x,'thbl:m[deg]' &
     &      ,1x,'phbl[deg]',1x,'pr[MeV/c]',1x,'isol')
  613 format(2x,'t1l[MeV]',3x,'th1l:m[deg]',1x,'ph1l[deg]',1x,'t2l[MeV]'   &
     &      ,3x,'th2l:m[MeV]',1x,'ph2l[deg]',1x,a2,'l',a7,1x,'thbl:m[deg]' &
     &      ,1x,'phbl[deg]',1x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  614 format(2x,'t1[MeV]',4x,'th1[deg]',4x,'ph1[deg]',2x,'t2[MeV]' &
     &      ,4x,'th2[deg]',4x,'ph2[deg]',2x,a2,a7,2x,'thb[deg]'    &
     &      ,4x,'phb[deg]',2x,'pr[MeV/c]',1x,'isol')
  615 format(2x,'t1[MeV]',4x,'th1[deg]',4x,'ph1[deg]',2x,'t2[MeV]' &
     &      ,4x,'th2[deg]',4x,'ph2[deg]',2x,a2,a7,2x,'thb[deg]'    &
     &      ,4x,'phb[deg]',2x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  616 format(2x,'t1[MeV]',4x,'th1:m[deg]',2x,'ph1[deg]',2x,'t2[MeV]' &
     &      ,4x,'th2:m[deg]',2x,'ph2[deg]',2x,a2,a7,2x,'thb:m[deg]'  &
     &      ,2x,'phb[deg]',2x,'pr[MeV/c]',1x,'isol')
  617 format(2x,'t1[MeV]',4x,'th1:m[deg]',2x,'ph1[deg]',2x,'t2[MeV]' &
     &      ,4x,'th2:m[deg]',2x,'ph2[deg]',2x,a2,a7,2x,'thb:m[deg]'  &
     &      ,2x,'phb[deg]',2x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  618 format(2x,'t1v[MeV]',3x,'th1v[deg]',3x,'ph1v[deg]',1x,'t2v[MeV]' &
     &      ,3x,'th2v[deg]',3x,'ph2v[deg]',1x,a2,'v',a7,1x,'thbv[deg]' &
     &      ,3x,'phbv[deg]',1x,'pr[MeV/c]',1x,'isol')
  619 format(2x,'t1v[MeV]',3x,'th1v[deg]',3x,'ph1v[deg]',1x,'t2v[MeV]' &
     &      ,3x,'th2v[deg]',3x,'ph2v[deg]',1x,a2,'v',a7,1x,'thbv[deg]' &
     &      ,3x,'phbv[deg]',1x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  620 format(2x,'t1v[MeV]',3x,'th1v:m[deg]',1x,'ph1v[deg]',1x,'t2v[MeV]'   &
     &      ,3x,'th2v:m[MeV]',1x,'ph2v[deg]',1x,a2,'v',a7,1x,'thbv:m[deg]' &
     &      ,1x,'phbv[deg]',1x,'pr[MeV/c]',1x,'isol')
  621 format(2x,'t1v[MeV]',3x,'th1v:m[deg]',1x,'ph1v[deg]',1x,'t2v[MeV]'   &
     &      ,3x,'th2v:m[MeV]',1x,'ph2v[deg]',1x,a2,'v',a7,1x,'thbv:m[deg]' &
     &      ,1x,'phbv[deg]',1x,'pr[MeV/c]',1x,'isol',1x,a20,2x,'Ay')

  701 format(1x,a2,'az'  ,a7,4x,'lgmd',a13)
  702 format(1x,a2,'az:m',a7,2x,'lgmd',a13)
  703 format(1x,a2,'ax'  ,a7,4x,'pxmd',a13)
  704 format(1x,a2,'ab'  ,a7,4x,'trmd',a13,2x,'trmd/',a2,a15)
  705 format(2x,a2,'a'   ,a7,4x,'tlmd',a13,2x,'tlmd/',a2,'^2',a15)

      return
end subroutine header
!-----------------------------------------------------------------------



!=======================================================================
subroutine mdcal(nkbabmax,dkbab,nkbazmax,fkbazmin,dkbaz,kunt,ddx)
      use consts
      use dims
      use kibmod
      implicit real*8 (a-h,o-z)
!=======================================================================
      real*8 ddx(nkzdim,nkbdim)
      real*8 fkbabms(nkbdim),fkbazms(nkbdim),trmd(nkbdim),pxmd(nkbdim)

!-----------------------------------------------------------------------

      if(kunt==0) then
        fkfac=1.0d0
       else if(kunt==1) then
        fkfac=hc
       else if(kunt==2) then
        fkfac=hc*1.0d-3
      end if

      if(kibpx+kibtr>0) then
        do ikbab=1,nkbabmax
         fkbab=(ikbab-1)*dkbab
         fkbabms(ikbab)=fkbab
         sum=0.0d0
         do ikbaz=1,nkbazmax
          sum=sum+ddx(ikbaz,ikbab)
         end do
         trmd(ikbab)=sum*dkbaz*fkfac
         if(kibtr>0)                                                     &
     &     write(kibtr,601) fkbab*fkfac,trmd(ikbab)*fkbab*fkfac*2.0d0*pi &
     &                     ,trmd(ikbab)*2.0d0*pi
  601     format(f12.5,2x,e13.5,6x,e13.5)
        end do

        if(kibpx>0) then
          fkbabmax=(nkbabmax-1)*dkbab
          do ikbax=1,nkbabmax
           fkbax=(ikbax-1)*dkbab
           sum=0.0d0
           do ikbay=1,2*nkbabmax-1
            fkbay=(ikbay-1)*dkbab-fkbabmax
            fkbab=sqrt(fkbax**2+fkbay**2)
            if(fkbab>fkbabmax) cycle
            comp=suphod(fkbab,fkbabms,dkbab,trmd,nkbabmax,nkbdim)
            if(comp>0.0d0) sum=sum+comp
           end do
           pxmd(ikbax)=sum*dkbab*fkfac
          end do
          do ikbax=1,2*nkbabmax-1
           fkbax=(ikbax-1)*dkbab-fkbabmax
           ikbaxc=nint(abs(fkbax)/dkbab)+1
           write(kibpx,602) fkbax*fkfac,pxmd(ikbaxc)
  602     format(f12.5,2x,e13.5)
          end do
        end if

      end if

      if(kibtl>0) then
        do ikbaz=1,nkbazmax
         fkbazms(ikbaz)=(ikbaz-1)*dkbaz+fkbazmin
        end do
        do ikbtl=1,nkbabmax
         fkbtl=(ikbtl-1)*dkbab
         sum=0.0d0
         do icx=-100,100
          cx=(icx-1)*1.0d-2
          if(icx==-100) cx=-1.0d0
          if(icx== 100) cx= 1.0d0
          sx=sqrt(1.0d0-cx**2)
          fkbaz=fkbtl*cx
          fkbab=fkbtl*sx
          if(fkbaz < fkbazms(1) .or. fkbaz > fkbazms(nkbazmax)) cycle
          comp=suphod2d(fkbaz,fkbazms,dkbaz,fkbab,fkbabms,dkbab &
     &                 ,ddx,nkbazmax,nkbabmax,nkzdim,nkbdim)
          if(comp>0.0d0) sum=sum+comp
         end do
         sum=sum*2.0d0*pi*1.0d-2
         write(kibtl,603) fkbtl*fkfac,sum*fkbtl**2*fkfac**2,sum
  603    format(f12.5,2x,e13.5,6x,e13.5)
        end do
      end if

      return

end subroutine mdcal
!=======================================================================
