      program nftincv10
c
c modified version
c
c hybrid, dual meshes, Finite volume edge based
c Navier_Stokes, laminar, incompressible
c NFT solver    using MPDATA and PK Smolarkiewicz's elliptic solver
c 2D solves for u,v,p (rho=const=1.)
c___________________________________________________________________
c Boundary condition flags
c * outer    1 - prescribes u & v to the free stream values
c * outer    7 - prescribes u (assuming it to be a no33rmal component)
c                to free strem value
c * symmetry 11 - normal velocity is forced to be zero
c * wall     21 - u=0 & v=0
c___________________________________________________________________
c
c author: Joanna Szmelter Loughborough University    Sept 2007
c
c___________________________________________________________________
c
c Created on 27/04/2018
c by Francesco Cocetta
c To check everything
c
c___________________________________________________________________
c
c NOT FOR DISTRIBUTION OR COMERCIAL USE WITHOUT AUTHOR'S PERMISSION
c
c___________________________________________________________________
cc      implicit double precision(a-h,o-z)
      open(5,file='mesh3d.d',access='sequential',status='old')
      open(35,file='mesh3d_c.d',access='sequential',status='old')
      open(45,file='mesh3d_cc.d',access='sequential',status='old')
      open(1,file='inputinc.d',access='sequential',status='old')
      open(8,file='input_mg_v3.d',access='sequential',status='old')
      open(12,file='inputre.d',access='sequential',status='old')
      open(55,file='match1.d',access='sequential',status='old')
      open(56,file='match2.d',access='sequential',status='old')
      open(11,file='r.d',access='sequential',status='unknown')
c
      open(7,file='trials.d',access='sequential',status='unknown')
      open(9,file='trials2.d',access='sequential',status='unknown')
      open(13,file='trials3.d',access='sequential',status='unknown')
c
      call cpu_time(cl_start)
c
      call input
      rewind 11
      call mpsolver
c
      call cpu_time(cl_end)
      print*, 'solution time is', cl_end-cl_start
c
      stop
      end
c
      subroutine mpsolver
c***********************************************************
c     NFT incompressible flow solver (usig MPDATA and conjugate
c     residua routines)
c***********************************************************
c
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      parameter (npoinx4=4*npoinx)
c
      common/per/iper
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/volume/vol(npoinx)
      common/itero/ itr,eps0,niter,nitsm,icount,eer,eem
      common/xyz/coord(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/solve/ntime,nrkstage
      common/sss/delt,adapt
      parameter (ielast=0)
      common/elast/ cfels,iflels
      dimension uuu(npoinx,3),uuu0(npoinx,3),uuu1(npoinx,3)
c    &                                      ,w(npoinx)
      dimension unk(npoinx,4),unk1(npoinx,4),dthedz(npoinx)
       data dthedz/npoinx*0./
      dimension rhs(npoinx),frc(npoinx,4),cour(npoinx)
      dimension prpot(npoinx),alpha(npoinx),unkdash(npoinx,3)
      common/coefficients/coefw(npoinx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
        data frc/npoinx4*0./
      data niter,nitsm,icount,eer,eem/3*0,2*0./
      common/ctherm/ rgas,cp,cap,strat,g,th00,
     &tt00,pr00,rh00,u00,v00,u0z,v0z
      DATA G/9.80616/
      common/pblfrc/ iblatt,hscl,vscl,hf00,prndt,cdrg
ccc  restart cccccc
      common/time/ntimet,nrstart,irestart,nre
c
      if(iper.eq.2)then
        call cyclicr(vol)
        call cyclicr2(vol)
        call cyclicr4(vol)
      endif
c
c read from inputflow.d      eps0=1.e-5
c read from inputflow.d      itr=100
c     f0=1.4584e-04 !coriolis parameter (=2xearth's omega)
      f0=1.4584e-03 !10x coriolis parameter (=2xearth's omega)
        xlat=20.
      iflels=ielast
      itraj=0
      implgw=1
      icorio=0
      ideep=1
      iblatt=0
      hscl=25.
      vscl=25.
      hf00=0.e-02
      prndt=0.42
      cdrg=0.0
c
      noutp=10
      nplot=200
c
      afree(1)=u0inf
      afree(2)=v0inf
      afree(3)=w0inf
      afree(4)=0.
      afree(5)=1.
c
      rh00=afree(5)
c      lipps=1
      lipps=0
      strat=1.e-05
      tetbar=300.
      tet00=tetbar
      tt00=tetbar
      th00=tetbar
c
      call rhprof
      call thprof
      call absorbers(alpha)
c     call noise(frc)

      pi=acos(-1.)
      tol=1.e-05

c set metric coeficients
      r=radious*float(isfera)+(1.-float(isfera))
      do ip=1,npoin
      y=coord(ip,2)
      cosa(ip)=cos(y)*float(isfera)+(1.-float(isfera))
      gmm(ip)=1.+coord(ip,3)/r*float(isfera)
      rh0(ip)=rh0(ip)*gmm(ip)**2*cosa(ip)
      coord(ip,1)=coord(ip,1)*r
      coord(ip,2)=coord(ip,2)*r
      enddo
c
c     compute time step beta and d metric factor
c     for the three meshes
      call beta_timestep_unstr(1)
      call beta_timestep_unstr(2)
      call beta_timestep_unstr(3)
      !
      !call compute_d()
      call compute_d_unstr(1)
      call compute_d_unstr(2)
      call compute_d_unstr(3)
      !call compute_d_c()
      !call compute_d_cc()
c
c set ambient wind (should be separate routine)
c     beta=-0.25*pi
      beta=0.
      cosb=cos(beta)
      sinb=sin(beta)
      do ip=1,npoin
      cosx=cos(coord(ip,1)/r)
      sinx=sin(coord(ip,1)/r)
      sina=sin(coord(ip,2)/r)
c------- construct geostrophic environmental state
cc    fcr3(i,j)=fcr0*(cobs*sina(i,j)-cosa(i,j)*cosx*sinb) !eulag example
      if(isfera.eq.1) then
      ue(ip)= afree(1)*(cosb*cosa(ip)+sina*cosx*sinb)
      ve(ip)=-afree(1)*sinx*sinb
      else
      ue(ip)= afree(1)
      ve(ip)= afree(2)
      endif

      if(lipps.eq.0)then
      the(ip)=tet00*(1.+strat*coord(ip,3))
      dthedz(ip) = tet00*strat*implgw                  ! pbl
      else
      the(ip)=th0(ip)
      dthedz(ip) = the(ip)*strat*implgw                ! pbl
      endif
      unk(ip,1)=ue(ip)
      unk(ip,2)=ve(ip)
      unk(ip,3)=0.2*frc(ip,1)
        unk(ip,4)=0.+0.001*frc(ip,1)
      unk1(ip,1)=unk(ip,1)/(gmm(ip)*cosa(ip))
      unk1(ip,2)=unk(ip,2)/gmm(ip)
      unk1(ip,3)=unk(ip,3)
      unk1(ip,4)=0.
      prpot(ip)=0.
      timestp(ip)=delt
      enddo
c numerical evaluation of dthedz
c     if(implgw.eq.1) then
c     call derivative(the,dthedz,3)
c     do ip=1,npoin
c     dthedz(ip)=dthedz(ip)/vol(ip)*implgw
c     enddo
c     endif

ccc      goto 1111
cc      call timecourant(uuu)
cc      call meshstats
c potential flow initialisation  -------------
      write(6,*)timestp(1)
      imp=0

      IF (nrstart .NE. 1)then
      call gcrk(prpot,unk1(1,1),unk1(1,2),unk1(1,3),
     &                 unk(1,1),unk(1,2),unk(1,3),imp)
      call prfr(prpot,unk1(1,1),unk1(1,2),unk1(1,3),
     &                 unk(1,1),unk(1,2),unk(1,3),imp)
      do ip=1,npoin
      uuu(ip,1)=unk(ip,1)
      uuu(ip,2)=unk(ip,2)
      uuu(ip,3)=unk(ip,3)
      unk(ip,1)=unk(ip,1)*gmm(ip)*cosa(ip)
      unk(ip,2)=unk(ip,2)*gmm(ip)
      enddo
      prsum=0.
      do ip=1,npoin
      prpot(ip)=-0.5*(unk(ip,1)*unk(ip,1)
     &+unk(ip,2)*unk(ip,2)+unk(ip,3)*unk(ip,3))/(0.5*timestp(ip))
      prsum=prsum+prpot(ip)
      enddo
      prsum=prsum/npoin
      do ip=1,npoin
      prpot(ip)=prpot(ip)-prsum
      prpot(ip)=0.
      frc(ip,1)=0.
      frc(ip,2)=0.
      frc(ip,3)=0.
      frc(ip,4)=0.
      gth0=g/th0(ip)
      coefw(ip)=1.+0.25*timestp(1)**2*gth0*dthedz(ip)
      coefw(ip)=1./coefw(ip)
      enddo
c----------------------------------------
      cfels=2./(300.*timestp(1))**2*iflels
c 1st order velocity prediction at n+0.5 in velprd
c end velprd
      if(itraj.eq.0) then
      do ip=1,npoin
      uuu(ip,1)=uuu(ip,1)*rh0(ip)
      uuu(ip,2)=uuu(ip,2)*rh0(ip)
      uuu(ip,3)=uuu(ip,3)*rh0(ip)
      unkdash(ip,1)=unk(ip,1)
      unkdash(ip,2)=unk(ip,2)
      unkdash(ip,3)=unk(ip,3)
      enddo
      else
      call velprd(unk,frc,uuu,uuu1,unkdash)
      endif
      call div(uuu(1,1),uuu(1,2),uuu(1,3),rhs,1)
      call courant(uuu,cour,1)
      call stats(unk)
      write(6,*)'end of initialisation, the time loop starts'


      else
      call Restart_In_old(uuu,uuu1,unkdash,unk,frc
     +   ,prpot,rhs,cour,coefw,nitsm,icount)


      write(*,*)'Restart from ',ntimet
      endif
c
c start iteration in time
      DO ITIME=1,NTIME
c     do ip=1,npoin
c     w(ip)=unk(ip,3)
c     enddo
ccccc
      do ip=1,npoin
       if(itraj.eq.0) then
      uuu0(ip,1)=uuu(ip,1)
      uuu0(ip,2)=uuu(ip,2)
      uuu0(ip,3)=uuu(ip,3)
       endif
      unk(ip,1)=unk(ip,1)+0.5*frc(ip,1)*timestp(ip)
      unk(ip,2)=unk(ip,2)+0.5*frc(ip,2)*timestp(ip)
      unk(ip,3)=unk(ip,3)+0.5*frc(ip,3)*timestp(ip)
      unk(ip,4)=unk(ip,4)+0.5*frc(ip,4)*timestp(ip)
        if(implgw.eq.0) unk(ip,4)=unk(ip,4)+the(ip)
      enddo
      if(isfera.eq.1)then
      call mpdatams(uuu,unk(1,1),1)
      call mpdatams(uuu,unk(1,2),2)
      call mpdatams(uuu,unk(1,3),3)
      call mpdatams(uuu,unk(1,4),4)
      endif
      if(isfera.eq.0.and.iper.eq.0)then
      call mpdatam(uuu,unk(1,1),1)
      call mpdatam(uuu,unk(1,2),2)
      call mpdatam(uuu,unk(1,3),3)
      call mpdatam(uuu,unk(1,4),4)
      endif
      if(iper.eq.2)then
      call mpdatamper(uuu,unk(1,1),1)
      call mpdatamper(uuu,unk(1,2),2)
      call mpdatamper(uuu,unk(1,3),3)
      call mpdatamper(uuu,unk(1,4),4)
      endif
c---------  -------------------------------------------------
      mtrimx=2
      do 292 mtri=1,mtrimx
c
      do 295 ip=1,npoin
       if(isfera.eq.1) then
      gmri=1./(gmm(ip)*r)
      tnga=sin(coord(ip,2)/r)/cosa(ip)
      frm1= gmri*tnga*(unkdash(ip,1)*unkdash(ip,2)-ue(ip)*ve(ip))
     &     -gmri*(unkdash(ip,1)*unkdash(ip,3))*ideep
      frm2=-gmri*tnga*(unkdash(ip,1)*unkdash(ip,1)-ue(ip)*ue(ip))
     &     -gmri*(unkdash(ip,2)*unkdash(ip,3))*ideep
      frm3= gmri*ideep*((unkdash(ip,1)*unkdash(ip,1)-ue(ip)*ue(ip))
     &     +(unkdash(ip,2)*unkdash(ip,2)-ve(ip)*ve(ip)) )
      fcr3=f0*sin(coord(ip,2)/r)
      fcr2=f0*cosa(ip)
       else
      frm1=0.
      frm2=0.
      frm3=0.
      fcr3=f0*sin(xlat)*icorio
      fcr2=f0*cos(xlat)*icorio
       endif
      if(mtri.eq.1.and.implgw.eq.0) unk(ip,4)=unk(ip,4)-the(ip)
      we=0.
      frm1=frm1+fcr3*(unkdash(ip,2)-ve(ip))
     &         -fcr2*(unkdash(ip,3)-we)*ideep
     &        -alpha(ip)*(unkdash(ip,1)-ue(ip))
      frm2=frm2-fcr3*(unkdash(ip,1)-ue(ip))
     &        -alpha(ip)*(unkdash(ip,2)-ve(ip))
      frm3=frm3+fcr2*(unkdash(ip,1)-ue(ip))*ideep
     &        -alpha(ip)*unkdash(ip,3)
     &+g*unk(ip,4)/th0(ip)                                  !ZMIANA
      frc(ip,1)=unk(ip,1)+0.5*timestp(ip)*frm1
      frc(ip,2)=unk(ip,2)+0.5*timestp(ip)*frm2
      frc(ip,3)=(unk(ip,3)+0.5*timestp(ip)*frm3)*coefw(ip)  !ZMIANA
      unk1(ip,1)=frc(ip,1)/(gmm(ip)*cosa(ip))
      unk1(ip,2)=frc(ip,2)/gmm(ip)
      unk1(ip,3)=frc(ip,3)
 295  continue
c
      imp=implgw
      call gcrk(prpot,unk1(1,1),unk1(1,2),unk1(1,3),
     &uuu(1,1),uuu(1,2),uuu(1,3),imp)
      call prfr(prpot,unk1(1,1),unk1(1,2),unk1(1,3),
     &uuu(1,1),uuu(1,2),uuu(1,3),imp)
c
      do 34 ip=1,npoin
      unkdash(ip,1)=uuu(ip,1)*gmm(ip)*cosa(ip)
      unkdash(ip,2)=uuu(ip,2)*gmm(ip)
      unkdash(ip,3)=uuu(ip,3)
      if(mtri.eq.mtrimx)then
        frc(ip,1)=(unkdash(ip,1)-frc(ip,1))*2./timestp(ip)
        frc(ip,2)=(unkdash(ip,2)-frc(ip,2))*2./timestp(ip)
        frc33=frc(ip,3)/coefw(ip)
     &       -(g*unk(ip,4)/th0(ip))*0.5*timestp(ip)
        frc(ip,3)=(unkdash(ip,3)-frc33)*2./timestp(ip)     !ZMIANA
      unk(ip,1)=unkdash(ip,1)
      unk(ip,2)=unkdash(ip,2)
      unk(ip,3)=unkdash(ip,3)
      frc(ip,4)=-unk(ip,3)*dthedz(ip)
      unk(ip,4)= unk(ip,4)+0.5*timestp(1)*frc(ip,4)
      endif
 34   continue
c
292   continue

      do 36 ip=1,npoin
       if(isfera.eq.1) then
      gmri=1./(gmm(ip)*r)
      tnga=sin(coord(ip,2)/r)/cosa(ip)
      frm1= gmri*tnga*(unk(ip,1)*unk(ip,2)-ue(ip)*ve(ip))
     &     -gmri*(unk(ip,1)*unk(ip,3))*ideep
      frm2=-gmri*tnga*(unk(ip,1)*unk(ip,1)-ue(ip)*ue(ip))
     &     -gmri*(unk(ip,2)*unk(ip,3))*ideep
      frm3= gmri*ideep*((unk(ip,1)*unk(ip,1)-ue(ip)*ue(ip))
     &     +(unk(ip,2)*unk(ip,2)-ve(ip)*ve(ip)) )
      fcr3=f0*sin(coord(ip,2)/r)
      fcr2=f0*cosa(ip)
        else
      frm1=0.
      frm2=0.
      frm3=0.
      fcr3=f0*sin(xlat)*icorio
      fcr2=f0*cos(xlat)*icorio
       endif
      frm1=frm1+fcr3*(unk(ip,2)-ve(ip))-fcr2*(unk(ip,3)-0.)*ideep
     &                                 -alpha(ip)*(unk(ip,1)-ue(ip))
      frm2=frm2-fcr3*(unk(ip,1)-ue(ip))-alpha(ip)*(unk(ip,2)-ve(ip))
      frm3=frm3+fcr2*(unk(ip,1)-ue(ip))*ideep-alpha(ip)*(unk(ip,3))
      frm4=                            -alpha(ip)*(unk(ip,4))
      frc(ip,1)=frc(ip,1)+frm1
      frc(ip,2)=frc(ip,2)+frm2
      frc(ip,3)=frc(ip,3)+frm3
      frc(ip,4)=frc(ip,4)+frm4

        if(iblatt.eq.1) then
         cfv=2-itraj
C  boundary layer forcing
        zatt=2.*exp(-coord(ip,3)/hscl)
        zatv=cfv*exp(-coord(ip,3)/vscl)
          frc(ip,4)=frc(ip,4)+(zatt/hscl)*hf00
          vnorm=cdrg*sqrt(unk(ip,1)**2+unk(ip,2)**2)*prndt
          frc(ip,1)=frc(ip,1)-(zatv/vscl)*vnorm*unk(ip,1)
          frc(ip,2)=frc(ip,2)-(zatv/vscl)*vnorm*unk(ip,2)
        endif

      if(itraj.eq.0) then
      uuu(ip,1)=1.5*uuu(ip,1)*rh0(ip)-0.5*uuu0(ip,1)
      uuu(ip,2)=1.5*uuu(ip,2)*rh0(ip)-0.5*uuu0(ip,2)
      uuu(ip,3)=1.5*uuu(ip,3)*rh0(ip)-0.5*uuu0(ip,3)
      rhoi=1./rh0(ip)
      uuu0(ip,1)=uuu0(ip,1)*rhoi*gmm(ip)*cosa(ip)
      uuu0(ip,2)=uuu0(ip,2)*rhoi*gmm(ip)
      uuu0(ip,3)=uuu0(ip,3)*rhoi
      unkdash(ip,1)=2.*unk(ip,1)-uuu0(ip,1)
      unkdash(ip,2)=2.*unk(ip,2)-uuu0(ip,2)
      unkdash(ip,3)=2.*unk(ip,3)-uuu0(ip,3)
      endif
  36  continue

      if(itraj.eq.1) then
       call velprd(unk,frc,uuu,uuu1,unkdash)
c second half after velprd
        if(iblatt.eq.1) then
      do ip=1,npoin
        zatv=exp(-coord(ip,3)/vscl)
          vnorm=cdrg*sqrt(unk(ip,1)**2+unk(ip,2)**2)*prndt
          frc(ip,1)=frc(ip,1)-(zatv/vscl)*vnorm*unk(ip,1)
          frc(ip,2)=frc(ip,2)-(zatv/vscl)*vnorm*unk(ip,2)
      enddo
        endif
      endif


cc       if(itime/noutp*noutp.eq.itime) then
      write(6,*)itime+ntimet, timestp(1),timestp(2)
      call div(uuu(1,1),uuu(1,2),uuu(1,3),rhs,1)
      call courant(uuu,cour,1)
      call stats(unk)
cc       endif
c      if(itime/nplot*nplot.eq.itime) call output(unk)
c
                   ENDDO
c     do ip=1,npoin
c     unk(ip,3)=0.5*(unk(ip,3)+w(ip))
c     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccc    Restart     ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      call Restart_Out(uuu,uuu1,unkdash,unk,frc,prpot,rhs,cour,
c     +                                               nitsm,icount,NTIME)
        call Restart_Out_old(uuu,uuu1,unkdash,unk,frc,prpot,rhs,cour,
     +                                       coefw,nitsm,icount,NTIME)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call output(unk)
      return
      end
c---------
c---------
      subroutine precon(rhs,p,r,iprc,elas,impfl)
cc      implicit double precision(a-h,o-z)
c***********************************************************
c preconditioner for the elliptic solver
c***********************************************************
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      common/xyz/coord(npoinx,3)
      common/xyz_c/coord_c(npoinx,3)
      common/edg/iedge(nedgex,2)
      common/edg_c/iedge_c(nedgex,2)
      !
      common/passelipt/ielast,epa,itrp,betap
      dimension rhs(npoinx),p(npoinx),r(npoinx),paux(npoinx)
      dimension rhs_c(npoinx),p_c(npoinx),r_c(npoinx)
      parameter (tol=1.e-08)
      !
      common/directions/nx,ny,nz
      common/directions_c/nx_c,ny_c,nz_c
      common/vcycle/mg_flag,omega,itrp1,itrp2,itrp_c1,itrp_c2,
     + itrp_cc
      common/smooth/sm_flag1,sm_flag2
      !
      dimension p_restr(npoinx)
      dimension p_prol(npoinx)
      ! for multigrid
      dimension rrhs(npoinx),rrhs_c(npoinx),rrhs_c_out(npoinx)
      dimension e(npoinx),e_c(npoinx)
      dimension pnew(npoinx),rnew(npoinx)
      !dimension d(npoinx)
      !
      common/subset/ie_subset(npoinx,2),sn_subset(npoinx,3),nsubset
      common/bound_point_flag/bou_flag(npoinx)
      !
      common/doperator/d1(npoinx),d2(npoinx),d3(npoinx)
      common/beta/betaunstr1,betaunstr2,betaunstr3
      ! vectors for trials
      dimension tr1(npoinx),tr2(npoinx),tr3(npoinx)
      !
      !
      ! coarse 2 grid
      common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
      common/edg_cc/iedge_cc(nedgex,2)
      common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      common/xyz_cc/coord_cc(npoinx,3)
      common/volume_cc/vol_cc(npoinx)
      common/wge_cc/sn_cc(nedgex,3)
      common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
      common/coin_cc/ncoin_cc,icoin_cc(nfacesx*2,2),irecog1_cc(nfacesx)
      common/terms_cc/sbb_cc(nfacesx,3)
      common/per2_cc/icoin2_cc(nfacesx,2),ncoin2_cc
      common/per4_cc/icoin4_cc(nfacesx,4),ncoin4_cc
      common/directions_cc/nx_cc,ny_cc,nz_cc
      !
      dimension rhs_cc(npoinx),p_cc(npoinx),r_cc(npoinx),rrhs_cc(npoinx)
      !
      dimension pstep(npoinx),pstep_c(npoinx),pstep_cc(npoinx)
c
c     this is not our case, iprc =1
c     in the mpsolver subroutine
      if(iprc.eq.0) then
      do ip=1,npoin
      p(ip)=rhs(ip)
      enddo
      return
      endif
c
c
c ---------------------------------------------------------
c
c     initialization of preconditioner
      do ip = 1,npoin
          p(ip)=0.
          r(ip)=0.
       enddo
c
c ----------------------------------------------------------
c
c    space for trials
c
*      do ip = 1,npoin_cc
*          p_c(ip)=0.
*          p_cc(ip)=1.
*      enddo
*      call prolongation_unstr(p_cc,p_c,2)
*      call prolongation_unstr(p_c,p,1)
*      do ic=1,npoin_c
*        write(9,*) ic, p_c(ic)
*        if (abs(p_c(ic)-1.).gt.1.e-8) then
*            print*, ic, "----check1 failed-----"
*        end if
*      end do
*      do ic=1,npoin
*        write(9,*) ic, p(ic)
*        if (abs(p(ic)-1.).gt.1.e-8) then
*            print*, ic, "----check2 failed-----"
*        end if
*      end do
*c
*      stop
c
c -----------------------------------------------------------
c
c     smoothing type for preconditioner
      if (sm_flag1.eq.1) then
          do ip=1,npoin
            pstep(ip)= (-1)/d1(ip)      ! Jacobi
          end do
      else if(sm_flag1.eq.2) then
          do ip=1,npoin
            pstep(ip)= betaunstr1       ! Richardson
          end do
      else if(sm_flag1.eq.3) then
          do ip=1,npoin
            pstep(ip)=betaunstr1/(1-d1(ip)*betaunstr1)        ! Mixed
          end do
      end if
c
c
c ---------------------------------------------------------
c
c   first smoothing on the finest grid
c
c
      do it = 1, itrp1
          call laplc(p,r,impfl)
          do ip = 1,npoin
              p(ip)=p(ip)+omega*pstep(ip)*(r(ip)-rhs(ip))
          end do
      end do
c
c ---------------------------------------------------------
c
c
c     MG enabled
      if (mg_flag.eq.1) then
c
c     Set up coarser grids
      do ip_c = 1,npoin_c
         p_c(ip_c)=0.  !0.
         r_c(ip_c)=0.  !0.
         e_c(ip_c)=0.
      enddo
      do ip_cc = 1,npoin_cc
         p_cc(ip_cc)=0.  !0.
         r_cc(ip_cc)=0.  !0.
      enddo
c
c     smoothing type for MG
      if (sm_flag1.eq.1) then
          do ic=1,npoin_c
            pstep_c(ic)= (-1)/d2(ic)      ! Jacobi
          end do
          do icc=1,npoin_cc
            pstep_cc(icc)= (-1)/d3(icc)  ! Jacobi
          end do
      else if(sm_flag1.eq.2) then
          do ic=1,npoin_c
            pstep_c(ic)= betaunstr2       ! Richardson
          end do
          do icc=1,npoin_cc
            pstep_cc(icc)= betaunstr3    ! Richardson
          end do
      else if(sm_flag1.eq.3) then
          do ic=1,npoin_c
            pstep_c(ic)=betaunstr2/(1-d2(ic)*betaunstr2)        ! Mixed
          end do
          do icc=1,npoin_cc
            pstep_cc(icc)=betaunstr3/(1-d3(icc)*betaunstr3)  ! Mixed
          end do
      end if
c
c
c -------------------------------------------------------------------
c
c residuals on the finest grid and restriction
      call laplc(p,r,impfl)
      do ip = 1, npoin
         rrhs(ip)=rhs(ip)-r(ip)
      end do
c
c      call restriction3D(rrhs,rrhs_c)
       call restriction_unstr(rrhs,rrhs_c,1)
c
c
c ---------------------------------------------------
c
      !
      ! relax on the coarse 1 grid
      !
      do ic = 1, npoin_c
          p_c(ic)= -(omega*pstep_c(ic)*rrhs_c(ic))
      end do
      do it=1, itrp_c1
          call laplc_c(p_c,r_c,impfl)
          do ic = 1, npoin_c
              p_c(ic)=p_c(ic)+omega*pstep_c(ic)*(r_c(ic)-rrhs_c(ic))
          enddo
      end do
      !
      ! ********************************************
      ! What I must supply to restriction in other stages?
      ! residuals on the coarse 1 grid?
      call laplc_c(p_c,r_c,impfl)
      do ip = 1, npoin_c
         rrhs_c_out(ip)=rrhs_c(ip)-r_c(ip)
      end do
      ! ********************************************
c
c -----------------------------------------------------------
c
      ! call the restriction on the coarse 2 grid
c      call restriction3D_cc(rrhs_c_out,rrhs_cc)
       call restriction_unstr(rrhs_c_out,rrhs_cc,2)
      !
      ! relax on the coarse 2 grid
      !
      do icc=1, npoin_cc
          p_cc(icc)= -(rrhs_cc(icc)*pstep_cc(icc)*omega)
      end do
      do  it=1, itrp_cc
          call laplc_cc(p_cc,r_cc,impfl)
          do icc = 1, npoin_cc
              p_cc(icc)=p_cc(icc)+omega*pstep_cc(icc)*
     + (r_cc(icc)-rrhs_cc(icc))
          enddo
      end do
c
c -----------------------------------------------------------
c
      ! prolongation to the coarse 1 grid and add the error
c      call prolongation3D_cc(p_cc,e_c)
      call prolongation_unstr(p_cc,e_c,2)
      do ic = 1,npoin_c
         p_c(ic)=p_c(ic)+e_c(ic)    !!!!
      enddo
      !
      !
      ! relaxation on the coarse 1 grid
      !
      do it=1, itrp_c2
            call laplc_c(p_c,r_c,impfl)
            do ic = 1, npoin_c
               p_c(ic)=p_c(ic)+omega*pstep_c(ic)*(r_c(ic)-rrhs_c(ic))
            enddo
      end do
      !
      ! prolongation to the finest grid
c      call prolongation3D(p_c,e)
      call prolongation_unstr(p_c,e,1)
      do ip = 1,npoin
         p(ip)=p(ip)+e(ip)  !!!!
      enddo
      !
      !
c ---------------------------------------------------------------------
      !
      ! final smoothing on the finest grid
      !
      do it = 1, itrp2  !0
          call laplc(p,r,impfl)
          do ip = 1, npoin
             p(ip)=p(ip)+omega*pstep(ip)*(r(ip)-rhs(ip))
          enddo
      end do
      end if ! flag on multigrid
c
c
c
      return
      end subroutine
c
c
c
      subroutine restriction_unstr(scalar,scalar_c,flg)
        !
        parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (npoinx4=4*npoinx)
        dimension scalar(npoinx),scalar_c(npoinx)
        !
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
        common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
        common/xyz/coord(npoinx,3)
        common/xyz_c/coord_c(npoinx,3)
        common/xyz_cc/coord_cc(npoinx,3)
        common/edg/iedge(nedgex,2)
        common/edg_c/iedge_c(nedgex,2)
        common/edg_cc/iedge_cc(nedgex,2)
        !
        common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
        common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
        common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
        common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
        common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
        common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
        !
        integer, dimension(npoinx) :: count_neigh,count_neigh_c,
     + count_neigh_cc
        integer, dimension(npoinx) :: match1,match2!,unmatch
        common/neighbours/count_neigh,neigh_list(npoinx,10)
        common/neighbours_c/count_neigh_c,neigh_list_c(npoinx,10)
        common/neighbours_cc/count_neigh_cc,neigh_list_cc(npoinx,10)
        common/match/match1,match2!,unmatch
        !
        ! flag to define the place of multigrid V-cyle I am
        ! and local variable definitions
        integer flg,n_coarse
        integer,dimension(npoinx) :: n_neigh_fine,match
        !
        ! if statement to define in what part of the V-cycle we are
        ! and what variable should be used
        if (flg.eq.1) then
            n_coarse = npoin_c
            n_neigh_fine = count_neigh
            match=match1
        else if (flg.eq.2) then
            n_coarse = npoin_cc
            n_neigh_fine = count_neigh_c
            match=match2
        end if
c
        scalar_c=0.
        ! loop on the coarse grid points
        do ic=1,n_coarse
            !
            ip=match(ic)!****map ic->ip
            !
            !
            ! for totally unstructured meshes
            do inei=1,n_neigh_fine(ip)
                !
                ! consider the mean value for every edge
                ! and later compute the global average
                scalar_c(ic)=scalar_c(ic)+(scalar(ip)+scalar(inei))/2
            end do
            scalar_c(ic)=scalar_c(ic)/n_neigh_fine(ip)
            !
            !
        end do

      end subroutine
c
c
c
      subroutine prolongation_unstr(scalar_c,scalar,flg)
        !
        parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (npoinx4=4*npoinx)
        dimension scalar(npoinx),scalar_c(npoinx),scalar_aux(npoinx)
        !
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
        common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
        common/xyz/coord(npoinx,3)
        common/xyz_c/coord_c(npoinx,3)
        common/xyz_cc/coord_cc(npoinx,3)
        common/edg/iedge(nedgex,2)
        common/edg_c/iedge_c(nedgex,2)
        common/edg_cc/iedge_cc(nedgex,2)
        !
        common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
        common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
        common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
        common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
        common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
        common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
        !
        integer, dimension(npoinx) :: count_neigh,count_neigh_c,
     + count_neigh_cc
        integer, dimension(npoinx) :: match1,match2!,unmatch
        common/neighbours/count_neigh,neigh_list(npoinx,10)
        common/neighbours_c/count_neigh_c,neigh_list_c(npoinx,10)
        common/neighbours_cc/count_neigh_cc,neigh_list_cc(npoinx,10)
        common/match/match1,match2!,unmatch
        !
        ! flag to define the place of multigrid V-cyle I am
        ! and local variable definitions
        integer flg,n_coarse
        integer,dimension(npoinx) :: n_neigh_fine,match
        integer flag_filling
        !
        ! if statement to define in what part of the V-cycle we are
        ! and what variable should be used
        if (flg.eq.1) then
            n_coarse = npoin_c
            npoin_fine = npoin
            n_neigh_fine = count_neigh
            match=match1
        else if (flg.eq.2) then
            n_coarse = npoin_cc
            npoin_fine = npoin_c
            n_neigh_fine = count_neigh_c
            match=match2
        end if
        !
        !
        nan=-99999
        scalar=0.
        scalar_aux=nan
        flag_filling=0
        !
        ! loop on the coarse grid points to transfer values of the
        ! coarse grid on the correspondent values of fine mesh,
        ! because it is an iterative process I need at least
        ! an auxiliary mesh
        !
        do ic=1,n_coarse
            !
            ip=match(ic)!****map ic->ip
            !
            scalar(ip)=scalar_c(ic)
            scalar_aux(ip)=scalar_c(ic)
        end do
        !
        ! for each point on the fine grid check for neighbours
        ! points and fill the empty spaces with the average value
        ! of the neighbours
  98  continue
        do ipa=1,npoin_fine
        !
c       ! with the other implementation the loop is only over
c       ! the points that not match with coarse mesh
c         do ipc=1,(npoin-npoin_c)
c            ipa=unmatch(ipc)
            if (scalar_aux(ipa).eq.nan) then
                !
                filled_neigh=0
                do ineigh=1,n_neigh_fine(ipa)
                    ! if statement to check if there are values
                    if (scalar_aux(ineigh).ne.nan) then
                      scalar(ipa)=scalar(ipa)+scalar(ineigh)
                      filled_neigh=filled_neigh+1
                    end if
                end do
                !
                ! check on division by zero and contemporary
                ! to not allocate all values,
                ! for Cartesian grid this do not matter
                !
                if (filled_neigh.ne.0) then
                    scalar(ipa)=scalar(ipa)/filled_neigh
                else
                    flag_filling=1
                end if
            end if
        end do
        !
        ! if the flag for unfilled values is 1, then put all
        ! the values on the auxiliary grid and restart the
        ! loop over primary grid points
        !
        if (flag_filling.eq.1) then
            print*, "All points are not filled"
            do ip=1,npoin_fine
              scalar_aux(ip)=scalar(ip)
            end do
            flag_filling=0
            go to 98
        end if
        !
        ! check if there still are NaN values
*        do ip=1,npoin
*           if(abs(scalar(ip)-0.).lt.1.e-03)then
*             print*, 'Point number', ip, 'is empty'
*           end if
*        end do
        !
      end subroutine
c
      subroutine beta_timestep_unstr(flg)
        !
        !
        ! This subroutine allows for the computation of the
        ! preconditioner timestep (beta) for fully unstructured meshes
        ! It should take into account for the different meshes
        !
        !
        parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (npoinx4=4*npoinx)
        common/domain/dom_x1,dom_x2,dom_y1,dom_y2,dom_z1,dom_z2
        !
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
        common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
        common/xyz/coord(npoinx,3)
        common/xyz_c/coord_c(npoinx,3)
        common/xyz_cc/coord_cc(npoinx,3)
        common/edg/iedge(nedgex,2)
        common/edg_c/iedge_c(nedgex,2)
        common/edg_cc/iedge_cc(nedgex,2)
        common/directions/nx,ny,nz
        !
        common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
        common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
        common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
        common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
        common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
        common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
        !
        common/terms/sbb(nfacesx,3)
        common/terms_c/sbb_c(nfacesx,3)
        common/terms_cc/sbb_cc(nfacesx,3)
        common/wge/sn(nedgex,3)
        common/wge_c/sn_c(nedgex,3)
        common/wge_cc/sn_cc(nedgex,3)
        common/volume/vol(npoinx)
        common/volume_c/vol_c(npoinx)
        common/volume_cc/vol_cc(npoinx)
        !
        integer, dimension(npoinx) :: count_neigh,count_neigh_c,
     + count_neigh_cc
        integer, dimension(npoinx) :: match1,match2!,unmatch
        common/neighbours/count_neigh,neigh_list(npoinx,10)
        common/neighbours_c/count_neigh_c,neigh_list_c(npoinx,10)
        common/neighbours_cc/count_neigh_cc,neigh_list_cc(npoinx,10)
        common/match/match1,match2!,unmatch
        !
        ! flag to define the place of multigrid V-cyle I am
        ! and local variable definitions
        integer flg,nedge_u,npoin_u
        integer,dimension(nedgex,2) :: iedge_u
        real,dimension(nedgex,3) :: sn_u
        real,dimension(npoinx) :: vol_u
        real,dimension(npoinx,3) :: coord_u
        !
        !output
        common/passelipt/ielast,epa,itrp,betap
        common/beta/betaunstr1,betaunstr2,betaunstr3
        !
        ! if statement to define in what part of the V-cycle we are
        ! and what variable should be used
        if (flg.eq.1) then
           nedge_u=nedge
           iedge_u=iedge
           sn_u=sn
           vol_u=vol
           coord_u=coord
           npoin_u=npoin
        else if (flg.eq.2) then
           nedge_u=nedge_c
           iedge_u=iedge_c
           sn_u=sn_c
           vol_u=vol_c
           coord_u=coord_c
           npoin_u=npoin_c
        else if (flg.eq.3) then
           nedge_u=nedge_cc
           iedge_u=iedge_cc
           sn_u=sn_cc
           vol_u=vol_cc
           coord_u=coord_cc
           npoin_u=npoin_cc
        end if
c
c
        !
        !
        ddx = 10.e+06
        find_min=0
        hill_height=1500.0
        !
        do ie=1,nedge_u
            !
            i1=iedge_u(ie,1)
            i2=iedge_u(ie,2)
            !
            ! IF statement to avoid points on the borders,
            ! because this points has got an associated value
            ! that is half of the internal points, making them
            ! negligible respect border points
            ! + last line to account for the presence of the hill
            ! = I consider the edge related to points above the hill
            !
            if ((coord_u(i1,1).ne.dom_x1).and.
     + (coord_u(i1,1).ne.dom_x2).and.
     + (coord_u(i1,2).ne.dom_y1).and.
     + (coord_u(i1,2).ne.dom_y2)
     + .and.(coord_u(i1,3).ne.dom_z1).and.
     + (coord_u(i1,3).ne.dom_z2).and.
     + (coord_u(i2,1).ne.dom_x1).and.
     + (coord_u(i2,1).ne.dom_x2)
     + .and.(coord_u(i2,2).ne.dom_y1).and.
     + (coord_u(i2,2).ne.dom_y2).and.
     + (coord_u(i2,3).ne.dom_z1)
     + .and.(coord_u(i2,3).ne.dom_z2).and.
     + (coord_u(i1,3).gt.hill_height).and.(coord_u(i2,3).gt.hill_height)
     + ) then
c
              !
              sx=sn_u(ie,1)
              sy=sn_u(ie,2)
              sz=sn_u(ie,3)
              !write(7,*) ie,sx,sy,sz,vol(i1),vol(i2)
              !
              hx=abs(sx)
              hy=abs(sy)
              hz=abs(sz)
              !
              h=max(hx,hy,hz)
c              if ((ddx.gt.vol_u(i1)/h).or.(ddx.gt.vol_u(i2)/h)) then
c                write(7,*) ie,i1,i2,sx,sy,sz,vol_u(i1),vol_u(i2),
c     + vol(i1)/h,vol(i2)/h
c              write(7,*)ie,(coord_u(i1,j),j=1,3),(coord_u(i2,j),j=1,3),h
c              write(7,*)i1, i2, vol_u(i1)/h,vol_u(i2)/h
c                find_min=ie
c              end if
              ddx=min(ddx,vol_u(i1)/h,vol_u(i2)/h)
            !
            end if
            !
        end do
        !
        !
        safety_factor=1 !4./3.
        beta_u=0.5*ddx*ddx*safety_factor
c        print*, beta_u, ddx
        !
        if (flg.eq.1) then
           betaunstr1=beta_u
        else if (flg.eq.2) then
           betaunstr2=beta_u
        else if (flg.eq.3) then
           betaunstr3=beta_u
        end if
c        stop
      end subroutine
c
c
      subroutine compute_d_unstr(flg)
        !
        parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (npoinx4=4*npoinx)
        !
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
        common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
        common/edg/iedge(nedgex,2)
        common/edg_c/iedge_c(nedgex,2)
        common/edg_cc/iedge_cc(nedgex,2)
        !
        common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
        common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
        common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
        common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
        common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
        common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
        !
        common/terms/sbb(nfacesx,3)
        common/terms_c/sbb_c(nfacesx,3)
        common/terms_cc/sbb_cc(nfacesx,3)
        common/wge/sn(nedgex,3)
        common/wge_c/sn_c(nedgex,3)
        common/wge_cc/sn_cc(nedgex,3)
        common/volume/vol(npoinx)
        common/volume_c/vol_c(npoinx)
        common/volume_cc/vol_cc(npoinx)
        !
        common/metric/cosa(npoinx),gmm(npoinx),radious
        common/coefficients/coefw(npoinx)
        common/tsteps/timestp(npoinx),istst,mrestart,time
        common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
        !
        ! local variables
        dimension dix(npoinx),diy(npoinx),diz(npoinx),d(npoinx)
        dimension ssumx(npoinx),ssumy(npoinx),ssumz(npoinx)
        !
        ! flag to define the place of multigrid V-cyle I am
        ! and local variable definitions
        integer flg,nedge_u,npoin_u,nbface_u
        integer,dimension(nedgex,2) :: iedge_u
        real,dimension(nedgex,3) :: sn_u
        real,dimension(npoinx) :: vol_u,rh0_u,rh0_c
        real,dimension(nfacesx,3) :: sb_u
        integer,dimension(nfacesx,1:2) :: kbface_u
        !
        !output
        common/doperator/d1(npoinx),d2(npoinx),d3(npoinx)
        !
        ! if statement to define in what part of the V-cycle we are
        ! and what variable should be used
        if (flg.eq.1) then
            npoin_u=npoin
            nedge_u=nedge
            iedge_u=iedge
            sn_u=sn
            vol_u=vol
            nbface_u=nbface
            kbface_u=kbface
            sb_u=sb
            rh0_u=rh0
        else if (flg.eq.2) then
            npoin_u=npoin_c
            nedge_u=nedge_c
            iedge_u=iedge_c
            sn_u=sn_c
            vol_u=vol_c
            nbface_u=nbface_c
            kbface_u=kbface_c
            sb_u=sb_c
            call restriction_unstr(rh0,rh0_u,1)
            ! note that if we deal with curvilinear domain
            ! we should restrict the geometrical factors
        else if (flg.eq.3) then
            npoin_u=npoin_cc
            nedge_u=nedge_cc
            iedge_u=iedge_cc
            sn_u=sn_cc
            vol_u=vol_cc
            nbface_u=nbface_cc
            kbface_u=kbface_cc
            sb_u=sb_cc
            call restriction_unstr(rh0,rh0_c,1)
            call restriction_unstr(rh0_c,rh0_u,2)
        end if
        !
        !
        ! Cycles to compute D_i;
        !
        do ip = 1,npoin_u
            dix(ip)=0.
            diy(ip)=0.
            diz(ip)=0.
            d(ip)=0.
            !
        end do
        !
        !
        ! Internal edges, including edges with one point on the faces
        do ied = 1,nedge_u
            !
            ! define the points related to each edge
            i1=iedge_u(ied,1)
            i2=iedge_u(ied,2)
            ! define the oriented surfaces
            sx=sn_u(ied,1)
            sy=sn_u(ied,2)
            sz=sn_u(ied,3)
            !write(7,*) ied, sx, sy, sz
            ! define the portion of oriented d for the points
            dix(i1)=dix(i1) + sx*sx*rh0_u(i2)/vol_u(i2)
            dix(i2)=dix(i2) + sx*sx*rh0_u(i1)/vol_u(i1)
            diy(i1)=diy(i1) + sy*sy*rh0_u(i2)/vol_u(i2)
            diy(i2)=diy(i2) + sy*sy*rh0_u(i1)/vol_u(i1)
            diz(i1)=diz(i1) + sz*sz*rh0_u(i2)/vol_u(i2)
            diz(i2)=diz(i2) + sz*sz*rh0_u(i1)/vol_u(i1)
            !
            !
        enddo
        !
        !
        ! Edges on the faces, it could be avoided in Cartesian mesh
        do ied = 1,nbface_u
            !
            ! define the points and surfaces related to each edge
            i1=kbface_u(ied,1)
            i2=kbface_u(ied,2)
            sx=sb_u(ied,1)
            sy=sb_u(ied,2)
            sz=sb_u(ied,3)
            ! define the portion of oriented d for the points
            ! rule is that i1->i2 is positive and i2->i1 is negative
            ! since there is a square, it does not matter
            dix(i1)=dix(i1) + sx*sx*rh0_u(i2)/vol_u(i2)
            dix(i2)=dix(i2) + sx*sx*rh0_u(i1)/vol_u(i1)
            diy(i1)=diy(i1) + sy*sy*rh0_u(i2)/vol_u(i2)
            diy(i2)=diy(i2) + sy*sy*rh0_u(i1)/vol_u(i1)
            diz(i1)=diz(i1) + sz*sz*rh0_u(i2)/vol_u(i2)
            diz(i2)=diz(i2) + sz*sz*rh0_u(i1)/vol_u(i1)
            !
            !
        enddo
        !
        !
        !impfl=0
        !
        do ip=1,npoin_u
            !
            ! following factors should be taken into account if
            ! we are dealing with
c            dix(ip)=dix(ip)/(gmm(ip)*cosa(ip))**2
c            diy(ip)=diy(ip)/(gmm(ip))**2
c            diz(ip)=diz(ip)*(coefw(ip)*impfl+(1-impfl))
            !
            d(ip)=-0.5*timestp(ip)*(dix(ip)+diy(ip)+diz(ip))/
     + (4*vol_u(ip)*rh0_u(ip))
            !
        end do
        !
        !
        ! reverse d into the right slot
        if (flg.eq.1) then
           d1=d
        else if (flg.eq.2) then
           d2=d
        else if (flg.eq.3) then
           d3=d
        end if
        !
        return
      end subroutine
c
c
c
      subroutine laplc(press,alap,impfl)
c***********************************************************
c computes laplacian for residuals and oses bc on residuals
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/coefficients/coefw(npoinx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume/vol(npoinx)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      dimension dpdxb(npoinx),dpdyb(npoinx),dpdzb(npoinx)
      dimension press(npoinx),alap(npoinx)
c
c
      do ip = 1,npoin
      alap(ip)=0.
      enddo
      call derivative(press,dpdxb,1)
      call derivative(press,dpdyb,2)
      call derivative(press,dpdzb,3)
      do ip=1,npoin
      dpdxb(ip)=timestp(ip)*dpdxb(ip)/vol(ip)
      dpdyb(ip)=timestp(ip)*dpdyb(ip)/vol(ip)
      dpdzb(ip)=timestp(ip)*dpdzb(ip)/vol(ip)
      dpdxb(ip)=0.5*dpdxb(ip)/(gmm(ip)*cosa(ip))**2
      dpdyb(ip)=0.5*dpdyb(ip)/gmm(ip)**2
      dpdzb(ip)=0.5*dpdzb(ip)*(coefw(ip)*impfl+(1-impfl))
      enddo
c
c note that the following "call bcv" clauds the transparency of
c codding.
c Here, it does not impse bc on velocity, but effectively impses
c bc conditions
c on the derivatives of residuals.
c 1) For boundaries for which all components
c of velocity are prescribed--- dpdx=0. and dpdy=0.;
c 2) For boundaries for which
c only one component of velocity is known ----dpdn=0.
c Other options for bc are possible but not implemented here.
      call bcv(dpdxb,dpdyb,dpdzb,0.)
c
c  loop over the internal edges
      do ip=1,npoin
      dpdxb(ip)=dpdxb(ip)*rh0(ip)
      dpdyb(ip)=dpdyb(ip)*rh0(ip)
      dpdzb(ip)=dpdzb(ip)*rh0(ip)
      enddo
      call derivative(dpdxb,alap,1)
      if(isfera.eq.1)then
      call derivativediv(dpdyb,dpdxb,2)
      else
      call derivative(dpdyb,dpdxb,2)
      endif
      call derivative(dpdzb,dpdyb,3)
c------------------------------------------
      do ip = 1,npoin
      alap(ip)=(alap(ip)+dpdxb(ip)+dpdyb(ip))/(vol(ip)*rh0(ip))
      enddo
c
      return
      end
c
c
      subroutine laplc_c(press,alap,impfl)
c***********************************************************
c computes laplacian for residuals and oses bc on residuals
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/coefficients/coefw(npoinx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume_c/vol_c(npoinx)
      common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      dimension dpdxb(npoinx),dpdyb(npoinx),dpdzb(npoinx)
      dimension press(npoinx),alap(npoinx)
      common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
      common/directions/nx,ny,nz
      common/directions_c/nx_c,ny_c,nz_c
      integer, dimension(npoinx)::match1,match2!,unmatch1,unmatch2
      common/match/match1,match2
c
c
      do ip = 1,npoin_c
      alap(ip)=0.
      enddo
      call derivative_c(press,dpdxb,1)
      call derivative_c(press,dpdyb,2)
      call derivative_c(press,dpdzb,3)
      do ip=1,npoin_c
      dpdxb(ip)=timestp(ip)*dpdxb(ip)/vol_c(ip)
      dpdyb(ip)=timestp(ip)*dpdyb(ip)/vol_c(ip)
      dpdzb(ip)=timestp(ip)*dpdzb(ip)/vol_c(ip)
      dpdxb(ip)=0.5*dpdxb(ip)/(gmm(ip)*cosa(ip))**2
      dpdyb(ip)=0.5*dpdyb(ip)/gmm(ip)**2
      dpdzb(ip)=0.5*dpdzb(ip)*(coefw(ip)*impfl+(1-impfl))
      enddo
c
c note that the following "call bcv" clauds the transparency of
c codding.
c Here, it does not impse bc on velocity, but effectively impses
c bc conditions
c on the derivatives of residuals.
c 1) For boundaries for which all components
c of velocity are prescribed--- dpdx=0. and dpdy=0.;
c 2) For boundaries for which
c only one component of velocity is known ----dpdn=0.
c Other options for bc are possible but not implemented here.
      call bcv_c(dpdxb,dpdyb,dpdzb,0.)
c
c
c     TRIALS
c      tol=1.e-7
c      write(7,*) "--------------"
c      do j=1,nboun_c
c        if((ibpoin_c(j,2).eq.7).or.(ibpoin_c(j,2).eq.8))then
c           ii=ibpoin_c(ib,1)
c           if (abs(dpdxb(ii)-0.).gt.tol) write(7,*)j, dpdxb(ii)
c        end if
c      enddo
c
c  loop over the internal edges
      do ip=1,npoin_c
        ! rh0 needs to be calibrated to account for the new grid.
        ! Here, for each point of the coarse grid is picked up
        ! the correspondent value of the finest grid
        !
        dpdxb(ip)=dpdxb(ip)*rh0(match1(ip))
        dpdyb(ip)=dpdyb(ip)*rh0(match1(ip))
        dpdzb(ip)=dpdzb(ip)*rh0(match1(ip))
      enddo
      call derivative_c(dpdxb,alap,1)
      if(isfera.eq.1)then
      call derivativediv(dpdyb,dpdxb,2)
      else
      call derivative_c(dpdyb,dpdxb,2)
      endif
      call derivative_c(dpdzb,dpdyb,3)
c------------------------------------------
      do ip = 1,npoin_c
        ! rh0 needs to be calibrated to account for the new grid.
        ! Here, for each point of the coarse grid is picked up
        ! the correspondent value of the finest grid
        !
        alap(ip)=(alap(ip)+dpdxb(ip)+dpdyb(ip))/
     + (vol_c(ip)*rh0(match1(ip)))
      enddo
c
      return
      end
c
c
c
      subroutine laplc_cc(press,alap,impfl)
c***********************************************************
c computes laplacian for residuals and oses bc on residuals
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/coefficients/coefw(npoinx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume_cc/vol_cc(npoinx)
      common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      dimension dpdxb(npoinx),dpdyb(npoinx),dpdzb(npoinx)
      dimension press(npoinx),alap(npoinx)
      common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
      integer, dimension(npoinx)::match1,match2!,unmatch1,unmatch2
      common/match/match1,match2
c
c
      do ip = 1,npoin_cc
      alap(ip)=0.
      enddo
      call derivative_cc(press,dpdxb,1)
      call derivative_cc(press,dpdyb,2)
      call derivative_cc(press,dpdzb,3)
      do ip=1,npoin_cc
      dpdxb(ip)=timestp(ip)*dpdxb(ip)/vol_cc(ip)
      dpdyb(ip)=timestp(ip)*dpdyb(ip)/vol_cc(ip)
      dpdzb(ip)=timestp(ip)*dpdzb(ip)/vol_cc(ip)
      dpdxb(ip)=0.5*dpdxb(ip)/(gmm(ip)*cosa(ip))**2
      dpdyb(ip)=0.5*dpdyb(ip)/gmm(ip)**2
      dpdzb(ip)=0.5*dpdzb(ip)*(coefw(ip)*impfl+(1-impfl))
      enddo
c
c note that the following "call bcv" clauds the transparency of
c codding.
c Here, it does not impse bc on velocity, but effectively impses
c bc conditions
c on the derivatives of residuals.
c 1) For boundaries for which all components
c of velocity are prescribed--- dpdx=0. and dpdy=0.;
c 2) For boundaries for which
c only one component of velocity is known ----dpdn=0.
c Other options for bc are possible but not implemented here.
      call bcv_cc(dpdxb,dpdyb,dpdzb,0.)
c
c
c     TRIALS
c      tol=1.e-7
c      write(7,*) "--------------"
c      do j=1,nboun_c
c        if((ibpoin_c(j,2).eq.7).or.(ibpoin_c(j,2).eq.8))then
c           ii=ibpoin_c(ib,1)
c           if (abs(dpdxb(ii)-0.).gt.tol) write(7,*)j, dpdxb(ii)
c        end if
c      enddo
c
c  loop over the internal edges
      do ip=1,npoin_cc
        ! rh0 needs to be calibrated to account for the new grid.
        ! Here, for each point of the coarse grid is picked up
        ! the correspondent value of the finest grid
        dpdxb(ip)=dpdxb(ip)*rh0(match1(match2(ip)))
        dpdyb(ip)=dpdyb(ip)*rh0(match1(match2(ip)))
        dpdzb(ip)=dpdzb(ip)*rh0(match1(match2(ip)))
      enddo
      call derivative_cc(dpdxb,alap,1)
      if(isfera.eq.1)then
      call derivativediv(dpdyb,dpdxb,2)
      else
      call derivative_cc(dpdyb,dpdxb,2)
      endif
      call derivative_cc(dpdzb,dpdyb,3)
c------------------------------------------
      do ip = 1,npoin_cc
        alap(ip)=(alap(ip)+dpdxb(ip)+dpdyb(ip))/
     + (vol_cc(ip)*rh0(match1(match2(ip))))
      enddo
c
      return
      end
c
c
      subroutine prfr(press,ustr,vstr,wstr,uout,vout,wout,impfl)
c***********************************************************
c computes pressure forces to determine velocities (uout,vout)
c at N+1 time level
c***********************************************************
c
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/coefficients/coefw(npoinx)
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume/vol(npoinx)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      dimension dpdxc(npoinx),dpdyc(npoinx),dpdzc(npoinx)
      dimension press(npoinx),ustr(npoinx),
     &vstr(npoinx),wstr(npoinx),uout(npoinx),vout(npoinx),wout(npoinx)
c
c
c  derivatives
      do ip = 1,npoin
      dpdxc(ip)=0.
      dpdyc(ip)=0.
      dpdzc(ip)=0.
      enddo
c
      call derivative(press,dpdxc,1)
      call derivative(press,dpdyc,2)
      call derivative(press,dpdzc,3)
      do ip=1,npoin
      dpdxc(ip)=dpdxc(ip)/(gmm(ip)*cosa(ip))**2/vol(ip)
      dpdyc(ip)=dpdyc(ip)/gmm(ip)**2/vol(ip)
      dpdzc(ip)=dpdzc(ip)/vol(ip)
      enddo
c
      do ip=1,npoin
      uout(ip)=ustr(ip)-0.5*timestp(ip)*dpdxc(ip)
      vout(ip)=vstr(ip)-0.5*timestp(ip)*dpdyc(ip)
      wout(ip)=wstr(ip)-0.5*timestp(ip)*dpdzc(ip)
     &        *(coefw(ip)*impfl+(1-impfl))
      enddo
      call bcv(uout,vout,wout,1.)
      return
      end
c-------------------
c-------------------
      subroutine bcv(uout,vout,wout,flag)
c***********************************************************
c applies bc for velocities (and residuals when called from
c Laplace)
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz/coord(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume/vol(npoinx)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      dimension uout(npoinx),vout(npoinx),wout(npoinx)
c  boundary conditions on velocities
c
      if(isfera.eq.1.or.iper.eq.2)then
      do ib=1,nboun
      ii=ibpoin(ib,1)
      if(ibpoin(ib,2).eq.11)then
      vn=uout(ii)*anp(ib,1)+vout(ii)*anp(ib,2)+wout(ii)*anp(ib,3)
      uout(ii)=uout(ii)-anp(ib,1)*vn
      vout(ii)=vout(ii)-anp(ib,2)*vn
      wout(ii)=wout(ii)-anp(ib,3)*vn
      endif
      enddo
      else
      do ib=1,nboun
      ii=ibpoin(ib,1)
      if(ibpoin(ib,2).eq.11)then
      vn=uout(ii)*anp(ib,1)+vout(ii)*anp(ib,2)+wout(ii)*anp(ib,3)
        uout(ii)=uout(ii)-anp(ib,1)*vn
        vout(ii)=vout(ii)-anp(ib,2)*vn
        wout(ii)=wout(ii)-anp(ib,3)*vn
        endif
      if(ibpoin(ib,2).eq.5)then
      vn=uout(ii)*anp(ib,1)+vout(ii)*anp(ib,2)+wout(ii)*anp(ib,3)
        uout(ii)=uout(ii)-anp(ib,1)*vn
        vout(ii)=vout(ii)-anp(ib,2)*vn
        wout(ii)=wout(ii)-anp(ib,3)*vn
        endif
      if(ibpoin(ib,2).eq.6)then
      vn=uout(ii)*anp(ib,1)+vout(ii)*anp(ib,2)+wout(ii)*anp(ib,3)
        uout(ii)=uout(ii)-anp(ib,1)*vn
        vout(ii)=vout(ii)-anp(ib,2)*vn
        wout(ii)=wout(ii)-anp(ib,3)*vn
        endif
      if(ibpoin(ib,2).eq.7)then
        uout(ii)=ue(ii)*flag
        endif
      if(ibpoin(ib,2).eq.8)then
        uout(ii)=ue(ii)*flag
        endif
      enddo
      endif
        return
        end
c--------------------------------------
c
      subroutine bcv_c(uout,vout,wout,flag)
c***********************************************************
c applies bc for velocities (and residuals when called from
c Laplace)
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz_c/coord_c(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume_c/vol_c(npoinx)
      common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
      common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      dimension uout(npoinx),vout(npoinx),wout(npoinx)
c  boundary conditions on velocities
c
      if(isfera.eq.1.or.iper.eq.2)then
      do ib=1,nboun_c
      ii=ibpoin_c(ib,1)
      if(ibpoin_c(ib,2).eq.11)then
      vn=uout(ii)*anp_c(ib,1)+vout(ii)*anp_c(ib,2)+wout(ii)*anp_c(ib,3)
      uout(ii)=uout(ii)-anp_c(ib,1)*vn
      vout(ii)=vout(ii)-anp_c(ib,2)*vn
      wout(ii)=wout(ii)-anp_c(ib,3)*vn
      endif
      enddo
      else
      do ib=1,nboun_c
      ii=ibpoin_c(ib,1)
      if(ibpoin_c(ib,2).eq.11)then
      vn=uout(ii)*anp_c(ib,1)+vout(ii)*anp_c(ib,2)+wout(ii)*anp_c(ib,3)
        uout(ii)=uout(ii)-anp_c(ib,1)*vn
        vout(ii)=vout(ii)-anp_c(ib,2)*vn
        wout(ii)=wout(ii)-anp_c(ib,3)*vn
        endif
      if(ibpoin_c(ib,2).eq.5)then
      vn=uout(ii)*anp_c(ib,1)+vout(ii)*anp_c(ib,2)+wout(ii)*anp_c(ib,3)
        uout(ii)=uout(ii)-anp_c(ib,1)*vn
        vout(ii)=vout(ii)-anp_c(ib,2)*vn
        wout(ii)=wout(ii)-anp_c(ib,3)*vn
        endif
      if(ibpoin_c(ib,2).eq.6)then
      vn=uout(ii)*anp_c(ib,1)+vout(ii)*anp_c(ib,2)+wout(ii)*anp_c(ib,3)
        uout(ii)=uout(ii)-anp_c(ib,1)*vn
        vout(ii)=vout(ii)-anp_c(ib,2)*vn
        wout(ii)=wout(ii)-anp_c(ib,3)*vn
        endif
      if(ibpoin_c(ib,2).eq.7)then
        uout(ii)=ue(ii)*flag
        endif
      if(ibpoin_c(ib,2).eq.8)then
        uout(ii)=ue(ii)*flag
        endif
      enddo
      endif
        return
        end
c
c
      subroutine bcv_cc(uout,vout,wout,flag)
c***********************************************************
c applies bc for velocities (and residuals when called from
c Laplace)
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz_cc/coord_cc(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume_cc/vol_cc(npoinx)
      common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
      common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      dimension uout(npoinx),vout(npoinx),wout(npoinx)
c  boundary conditions on velocities
c
      if(isfera.eq.1.or.iper.eq.2)then
      do ib=1,nboun_cc
      ii=ibpoin_cc(ib,1)
      if(ibpoin_cc(ib,2).eq.11)then
      vn=uout(ii)*anp_cc(ib,1)+vout(ii)*anp_cc(ib,2)+
     + wout(ii)*anp_cc(ib,3)
      uout(ii)=uout(ii)-anp_cc(ib,1)*vn
      vout(ii)=vout(ii)-anp_cc(ib,2)*vn
      wout(ii)=wout(ii)-anp_cc(ib,3)*vn
      endif
      enddo
      else
      do ib=1,nboun_cc
      ii=ibpoin_cc(ib,1)
      if(ibpoin_cc(ib,2).eq.11)then
      vn=uout(ii)*anp_cc(ib,1)+vout(ii)*anp_cc(ib,2)+
     + wout(ii)*anp_cc(ib,3)
        uout(ii)=uout(ii)-anp_cc(ib,1)*vn
        vout(ii)=vout(ii)-anp_cc(ib,2)*vn
        wout(ii)=wout(ii)-anp_cc(ib,3)*vn
        endif
      if(ibpoin_cc(ib,2).eq.5)then
      vn=uout(ii)*anp_cc(ib,1)+vout(ii)*anp_cc(ib,2)+
     + wout(ii)*anp_cc(ib,3)
        uout(ii)=uout(ii)-anp_cc(ib,1)*vn
        vout(ii)=vout(ii)-anp_cc(ib,2)*vn
        wout(ii)=wout(ii)-anp_cc(ib,3)*vn
        endif
      if(ibpoin_cc(ib,2).eq.6)then
      vn=uout(ii)*anp_cc(ib,1)+vout(ii)*anp_cc(ib,2)+
     + wout(ii)*anp_cc(ib,3)
        uout(ii)=uout(ii)-anp_cc(ib,1)*vn
        vout(ii)=vout(ii)-anp_cc(ib,2)*vn
        wout(ii)=wout(ii)-anp_cc(ib,3)*vn
        endif
      if(ibpoin_cc(ib,2).eq.7)then
        uout(ii)=ue(ii)*flag
        endif
      if(ibpoin_cc(ib,2).eq.8)then
        uout(ii)=ue(ii)*flag
        endif
      enddo
      endif
        return
        end
c
c
      subroutine gcrk(p,u1,v1,w1,u2,v2,w2,impfl)
c***********************************************************
c Piotr Smolarkiewicz elliptic solver
c***********************************************************
      common/passelipt/ielast,epa,itrp,betap
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      parameter (lordx=9,nn=npoinx,lord=9)
c      parameter (nplt=100)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/itero/ itr,eps0,niter,nitsm,icount,eer,eem
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/elast/ cfels,iflels
      dimension p(nn),u1(nn),v1(nn),u2(nn),v2(nn),w1(nn),w2(nn),
     &r(nn),qr(nn),ar(nn),
     *x(nn,lordx),ax(nn,lordx),ax2(lordx),axar(lordx),del(lordx)
      dimension p0(npoinx)
c      dimension err(0:nplt)
convergence test modes *********************************************
      logical ctest
      data ctest/.false./
cc    data ctest/.true./
c      nplt=100
c      if(ctest) then
c      itr= 100/lord
c      ner= 1
c      snorm=1./float(nn)
c      endif
convergence test modes *********************************************
cc in inputflow.d now      eps0=1.e-15
      eps=eps0/timestp(1)
      nml=npoin
      epa=1.e-38
      nlc=0
      itmn=3
c      iprc=0
      iprc=1
      if(iflels.eq.1) then
       do k=1,nml
        p0(k)=p(k)
       enddo
      endif

      do k=1,nml
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do i=1,lord
       do k=1,nml
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call prfr(p,u1,v1,w1,u2,v2,w2,impfl)
      call div(u2,v2,w2,r,-1)
        if(iflels.eq.1) then
         do k=1,nml
          r(k)=r(k)-elas*(p(k)-p0(k))
         enddo
        endif
        call precon(r,qr,ar,iprc,elas,impfl)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nml
      eer0=eer0+qr(k)**2
      eem0=amax1(eem0,abs(qr(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=amax1(eer0,epa)
      eem0=amax1(eem0,epa)
      rl20=amax1(rl20,epa)
convergence test modes **********************************************
c      if(ctest) then
c      do ier=0,nplt
c      err(ier)=eps
c      enddo
c      eer=-1.e15
c      do 3 k=1,nml
c    3 eer=amax1(eer,abs(r(k)))
c      err(0)=eer
c      print 300,  err(0)
c  300 format(4x,e11.4,' residual error at it=1')
c      endif
convergence test modes **********************************************
       do k=1,nml
        x(k,1)=qr(k)
       enddo
      call laplc(x(1,1),ax(1,1),impfl)
        if(iflels.eq.1) then
         do k=1,nml
          ax(k,1)=ax(k,1)-elas*x(k,1)
         enddo
        endif
c     if(joasia.eq.0)then
c      itr1=itr
c     itr=itr*6
c      joasia=joasia+1
c     else
c      itr=itr1
c     endif
      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nml
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=amax1(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
         do k=1,nml
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=amax1(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
         enddo
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=amax1(rl2,epa)
       call precon(r,qr,ar,iprc,elas,impfl)
       call laplc(qr,ar,impfl)
        if(iflels.eq.1) then
         do k=1,nml
          ar(k)=ar(k)-elas*r(k)
         enddo
        endif
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nml
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
c         del(ii)=amax1(del(ii),0.5)
         enddo
        if(i.lt.lord) then
          do k=1,nml
            x(k,i+1)=qr(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nml
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nml
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nml
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
convergence test modes **********************************************
c      if(ctest) then
c      if(nlc/ner*ner.eq.nlc) then
c      ier=nlc/ner
c      eer=-1.e15
c      do 50 k=1,nml
c   50 eer=amax1(eer,abs(r(k)))
c      err(ier)=eer
c      endif
c      endif
convergence test modes **********************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nml
      eer=eer+qr(k)**2
      eem=amax1(eem,abs(qr(k)))
      enddo
      eer=eer/eer0
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

convergence test modes **********************************************
c      if(ctest) then
c      print 301, (err(ier),ier=1,nplt,1)
c  301 format(4x,5e11.4)
c      endif
convergence test modes **********************************************
      return
      end
c
      subroutine courant(uuu,cour,iflag)
c***********************************************************
c provides diagnostic information
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension uuu(npoinx,3),cour(npoinx)
c
      do ip = 1,npoin
      cour(ip)=0.
      enddo
c  loop over the internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
      rhav=0.5*(rh0(i1)+rh0(i2))
      aun=(uuu(i1,1)+uuu(i2,1))*0.5*sx
     &   +(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &   +(uuu(i1,3)+uuu(i2,3))*0.5*sz
      aun=aun/rhav
      if(iflag.eq.1) then
       cour(i1)=cour(i1)+abs(aun)
       cour(i2)=cour(i2)+abs(aun)
      endif
      if(iflag.eq.2) then
       cour(i1)=cour(i1)+aun
       cour(i2)=cour(i2)-aun
      endif
      if(iflag.eq.3) then
       cour(i1)=cour(i1)+amax1(0.,aun)
       cour(i2)=cour(i2)-amin1(0.,aun)
      endif
      enddo
c
c  loop over boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
      rhav=0.5*(rh0(i1)+rh0(i2))
      aun=(uuu(i1,1)+uuu(i2,1))*0.5*sx
     &   +(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &   +(uuu(i1,3)+uuu(i2,3))*0.5*sz
      aun=aun/rhav
      if(iflag.eq.1) then
       cour(i1)=cour(i1)+abs(aun)
       cour(i2)=cour(i2)+abs(aun)
      endif
      if(iflag.eq.2) then
       cour(i1)=cour(i1)+aun
       cour(i2)=cour(i2)-aun
      endif
      if(iflag.eq.3) then
       cour(i1)=cour(i1)+amax1(0.,aun)
       cour(i2)=cour(i2)-amin1(0.,aun)
      endif
      enddo
c  loop over boundary faces edges
      if(iper.eq.2)goto 987
      do  ied = 1,nboun
      i1=ibpoin(ied,1)
      sx=sbb(ied,1)
      sy=sbb(ied,2)
      sz=sbb(ied,3)
      rhav=rh0(i1)
      aun=uuu(i1,1)*sx
     &   +uuu(i1,2)*sy
     &   +uuu(i1,3)*sz
      aun=aun/rhav
      if(iflag.eq.1) then
       cour(i1)=cour(i1)+abs(aun)
      endif
      if(iflag.eq.2) then
       cour(i1)=cour(i1)+aun
      endif
      if(iflag.eq.3) then
       cour(i1)=cour(i1)+amax1(0.,aun)
      endif
      enddo
 987  continue
      do ip = 1,npoin
        cour(ip)=cour(ip)*timestp(ip)/vol(ip)/rh0(ip)
      enddo
c print diagnostics
      crnmx=-1.e9
      crnmn= 1.e9
      crnav=0.
        crnsd=0.
      do ip=1,npoin
        crnav=crnav+cour(ip)
        crnmx=amax1(crnmx,cour(ip))
        crnmn=amin1(crnmn,cour(ip))
       if(crnmx.eq.cour(ip))ic=ip
      enddo
        crnav=crnav/npoin
        do ip=1,npoin
        crnsd=crnsd+(cour(ip)-crnav)**2
        enddo
        crnsd=sqrt(crnsd/npoin)
      write(6,*)ic
      write(6,205)crnmx,crnmn,crnav,crnsd
  205 format(1x,'crnmx,crnmn,crnav,crnsd:',4e11.4)
      return
      end
c--------------------------------------
      subroutine div(ua,va,wa,rhs,iflg)
c***********************************************************
c computes divergence (and supplies diagnostics)
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz/coord(npoinx,3)
      common/itero/ itr,eps0,niter,nitsm,icount,eer,eem
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension dudx(npoinx),dvdy(npoinx),dwdz(npoinx)
      dimension ua(npoinx),va(npoinx),wa(npoinx)
      dimension rhs(npoinx),ub(npoinx),vb(npoinx),wb(npoinx)
c
c
      if(iflg.eq.-1) then
      do ip=1,npoin
      ub(ip)=ua(ip)*rh0(ip)
      vb(ip)=va(ip)*rh0(ip)
      wb(ip)=wa(ip)*rh0(ip)
      enddo
      else
      do ip=1,npoin
      ub(ip)=ua(ip)
      vb(ip)=va(ip)
      wb(ip)=wa(ip)
      enddo
      endif
c  derivatives
      call derivative(ub,dudx,1)
      if(isfera.eq.1)then
      call derivativediv(vb,dvdy,2)
      else
      call derivative(vb,dvdy,2)
      endif
      call derivative(wb,dwdz,3)
c------------------------------------------
      do ip = 1,npoin
      rhs(ip)=iflg*(dudx(ip)/vol(ip)+dvdy(ip)/vol(ip)
     *             +dwdz(ip)/vol(ip))/rh0(ip)
      enddo
c
      if(iflg.eq.1) then
      divmx=-1.e15
      divmn= 1.e15
      divav=0.
      do ip =1,npoin
      divmx=amax1(divmx,rhs(ip)*timestp(ip))
      divmn=amin1(divmn,rhs(ip)*timestp(ip))
      divav=divav+rhs(ip)*timestp(ip)
      enddo
      divav=divav/float(npoin)
      divsd=0.
      do ip=1,npoin
      divsd=divsd+(rhs(ip)*timestp(ip)-divav)**2
      enddo
      divsd=sqrt(divsd/npoin)
      nitav=nitsm/max0(icount,1)
      print 205, divmx,divmn,divav,divsd,
     +eer,eem,niter,nitav
  205 format(1x,'dvmx,dvmn,dvav,dvsd:',4e11.4/
     .1x,' eer,eem:',2e11.4,1x,'niter,nitav:',4i4)
      endif
      return
      end
c---------------------------------------------
      subroutine input
c***********************************************************
c     reads geometrical and flow  data
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/domain/dom_x1,dom_x2,dom_y1,dom_y2,dom_z1,dom_z2
c
      common/per/iper
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/per2/  icoin2(nfacesx,2),ncoin2
      common/per4/  icoin4(nfacesx,4),ncoin4
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/terms/sbb(nfacesx,3)
      common/wge/sn(nedgex,3)
      common/itero/ itr,eps0,niter,nitsm,icount,eer,eem
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/sss/delt,adapt
      common/ss/switch(npoinx),relax,switch1(npoinx)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      common/solve/ntime,nrkstage
      common/cflpass/cfl
      common/absspec/dxabl,dxabr,dyabd,dyabu,dzabd,dzabu,
     &towxl,towxr,towyd,towyu,towzd,towzu,irl
      common/passelipt/ielast,epa,itrp,betap
      common/volume/vol(npoinx)
      common/sf/isfera
ccc  restart cccccc
      common/time/ntimet,nrstart,irestart,nre
      !
      ! common variables for reading the coarse mesh
      !
      common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
      common/edg_c/iedge_c(nedgex,2)
      common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      common/xyz_c/coord_c(npoinx,3)
      common/volume_c/vol_c(npoinx)
      common/wge_c/sn_c(nedgex,3)
      common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
      common/coin_c/ncoin_c,icoin_c(nfacesx*2,2),irecog1_c(nfacesx)
      common/terms_c/sbb_c(nfacesx,3)
      common/per2_c/icoin2_c(nfacesx,2),ncoin2_c
      common/per4_c/icoin4_c(nfacesx,4),ncoin4_c
      !
      ! variables for the second coarse mesh
      !
      common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
      common/edg_cc/iedge_cc(nedgex,2)
      common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      common/xyz_cc/coord_cc(npoinx,3)
      common/volume_cc/vol_cc(npoinx)
      common/wge_cc/sn_cc(nedgex,3)
      common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
      common/coin_cc/ncoin_cc,icoin_cc(nfacesx*2,2),irecog1_cc(nfacesx)
      common/terms_cc/sbb_cc(nfacesx,3)
      common/per2_cc/icoin2_cc(nfacesx,2),ncoin2_cc
      common/per4_cc/icoin4_cc(nfacesx,4),ncoin4_cc
      !
      ! -------------
      ! Following variables are defined for unstructured
      ! subroutines, storing matching and neighbours values
      integer, dimension(npoinx)::count_neigh,count_neigh_c,
     + count_neigh_cc
      integer, dimension(npoinx)::match1,match2!,unmatch1,unmatch2
      common/neighbours/count_neigh,neigh_list(npoinx,10)
      common/neighbours_c/count_neigh_c,neigh_list_c(npoinx,10)
      common/neighbours_cc/count_neigh_cc,neigh_list_cc(npoinx,10)
      common/match/match1,match2!,unmatch1,unmatch2
      !
      ! Multigrid cycle characteristics
      common/vcycle/mg_flag,omega,itrp1,itrp2,itrp_c1,itrp_c2,
     + itrp_cc
      common/smooth/sm_flag1,sm_flag2
      integer match_flag
c
c
c manually inserting the boundaries of the domain
c
      dom_x1=0.00000
      dom_x2=15000.00000
      dom_y1=-6000.00000
      dom_y2=+6000.00000
      dom_z1=0.00000
      dom_z2=6000.00000
c
c
c read file for number of Richardson iterations in the fine and coarse (2) grid
c
      read(8,*)
      read(8,*) mg_flag
      read(8,*)
      read(8,*)
      read(8,*) sm_flag1, sm_flag2
      read(8,*)
      read(8,*) omega
      read(8,*)
      read(8,*)
      read(8,*) itrp1, itrp2, itrp_c1, itrp_c2, itrp_cc
      read(8,*)
      read(8,*) match_flag
c
c
c
c
c read a control file
c read flow parameters
      read(1,*)
c ntime = number of time iterations
c nrkstage= iord= number of MPDATA passes
c
c istst=0 for global time step
c istst=1 for local timestep
c istst=2 const delt prescribed in the inputflow.d
c
c to restart set mrestart=1 and change the name from results.d
cobtained in the
c previous run to results0.d
c
      read(1,*)ntime,nrkstage,istst,mrestart,isfera,iper
      read(1,*)
      read(1,*)cfl,delt,adapt,relax
      read(1,*)
c read free stream pressure, velocities and Re number
      read(1,*)p0inf,u0inf,v0inf,w0inf,Re              ! liu:    ,w0inf  added
      write(6,'(a)')' p0inf,u0inf,v0inf,w0inf,Re'
      write(6,*)p0inf,u0inf,v0inf,w0inf
c
c      write(12,'(a)')' p0inf,u0inf,v0inf,w0inf,Re'   ! liu
c      write(12,*)p0inf,u0inf,v0inf,w0inf
c
c read specification for elliptic solver
      read(1,*)
      read(1,*)
      read(1,*)ielast,eps0,itr,lord
      read(1,*)
      read(1,*)itrp,betap
c read specifications for absorbers
      read(1,*)
      read(1,*)irl
      if(irl.eq.1)then
      read(1,*)
      read(1,*)dxabl,dxabr,dyabd,dyabu,dzabd,dzabu
      read(1,*)
      read(1,*)towxl,towxr,towyd,towyu,towzd,towzu
      endif
c
c restart file
      READ(12,*)
      READ(12,*)nrstart
      READ(12,*)
      READ(12,*)ntimet
      READ(12,*)
      READ(12,*)irestart,nre
c
c     radious=6371.22e+03
c     radious=30.e+03
      radious=63.6620e+03
      if(isfera.eq.0)radious=1.
      r=radious
c
c Definition of the number of points for each mesh
c
c
c read mesh first
c Finest mesh
c
      read(5,*)npoin, nedge, nbface, nbpoin, nboun
      write(6,*)'Primary mesh',npoin,nedge,nbface,nbpoin,nboun
      if(npoinx .lt. npoin)then             ! liu:  check maximum
        write(*,*)'npoinx < npoin'
        write(*,'(" npoinx =",i8)') npoinx
        write(*,'(" npoin =",i8)') npoin
        stop
      endif
      if(nedgex .lt. nedge)then
        write(*,*)'nedgex < nedge'
        write(*,'(" nedgex =",i8)') nedgex
        write(*,'(" nedge =",i8)') nedge
        stop
      endif
      if(nfacesx .lt. nbface)then
        write(*,*)'nfacesx < nbface'
        write(*,'(" nfacesx =",i8)') nfacesx
        write(*,'(" nbface =",i8)') nbface
        stop
      endif
c
c Coarse 1 mesh
c
      read(35,*)npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      write(6,*)'Coarse 1 mesh',npoin_c,nedge_c,
     + nbface_c,nbpoin_c,nboun_c
      if(npoinx .lt. npoin_c)then             ! liu:  check maximum
        write(*,*)'npoinx < npoin'
        write(*,'(" npoinx =",i8)') npoinx
        write(*,'(" npoin =",i8)') npoin_c
        stop
      endif
      if(nedgex .lt. nedge_c)then
        write(*,*)'nedgex < nedge'
        write(*,'(" nedgex =",i8)') nedgex
        write(*,'(" nedge =",i8)') nedge_c
        stop
      endif
      if(nfacesx .lt. nbface_c)then
        write(*,*)'nfacesx < nbface'
        write(*,'(" nfacesx =",i8)') nfacesx
        write(*,'(" nbface =",i8)') nbface_c
        stop
      endif
c
c Coarse 2 mesh
c
      read(45,*)npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      write(6,*)'Coarse 2 mesh',npoin_cc,nedge_cc,
     + nbface_cc,nbpoin_cc,nboun_cc
      if(npoinx .lt. npoin_cc)then             ! liu:  check maximum
        write(*,*)'npoinx < npoin'
        write(*,'(" npoinx =",i8)') npoinx
        write(*,'(" npoin =",i8)') npoin_cc
        stop
      endif
      if(nedgex .lt. nedge_cc)then
        write(*,*)'nedgex < nedge'
        write(*,'(" nedgex =",i8)') nedgex
        write(*,'(" nedge =",i8)') nedge_cc
        stop
      endif
      if(nfacesx .lt. nbface_cc)then
        write(*,*)'nfacesx < nbface'
        write(*,'(" nfacesx =",i8)') nfacesx
        write(*,'(" nbface =",i8)') nbface_cc
        stop
      endif
c
c read coordinates
c finest mesh
c
      do ip=1,npoin
        read(5,*)ii,(coord(ip,i),i=1,3),vol(ip)
      vol(ip)=vol(ip)*r**2
      enddo
c
c coarse 1 mesh
c
      do ip_c=1,npoin_c
        read(35,*)ii,(coord_c(ip_c,i),i=1,3),vol_c(ip_c)
      vol_c(ip_c)=vol_c(ip_c)*r**2
      enddo
c
c coarse 2 mesh
c
      do ip_cc=1,npoin_cc
        read(45,*)ii,(coord_cc(ip_cc,i),i=1,3),vol_cc(ip_cc)
      vol_cc(ip_cc)=vol_cc(ip_cc)*r**2
      enddo
c
c read edges
c finest mesh
c
      do ie=1,nedge
        read(5,*)iie,(iedge(ie,i),i=1,2)
      enddo
      do ie=1,nedge
        read(5,*)iie, (sn(ie,i),i=1,3)
      sn(ie,1)=sn(ie,1)*r
      sn(ie,2)=sn(ie,2)*r
      sn(ie,3)=sn(ie,3)*r*r
      enddo
c
c coarse 1 mesh
c
      do ie_c=1,nedge_c
        read(35,*)iie_c, (iedge_c(ie_c,i),i=1,2)
      enddo
      do ie_c=1,nedge_c
        read(35,*)iie_c, (sn_c(ie_c,i),i=1,3)
      sn_c(ie_c,1)=sn_c(ie_c,1)*r
      sn_c(ie_c,2)=sn_c(ie_c,2)*r
      sn_c(ie_c,3)=sn_c(ie_c,3)*r*r
      enddo
c
c coarse 2 mesh
c
      do ie_cc=1,nedge_cc
        read(45,*)iie_cc, (iedge_cc(ie_cc,i),i=1,2)
      enddo
      do ie_cc=1,nedge_cc
        read(45,*)iie_cc, (sn_cc(ie_cc,i),i=1,3)
      sn_cc(ie_cc,1)=sn_cc(ie_cc,1)*r
      sn_cc(ie_cc,2)=sn_c(ie_cc,2)*r
      sn_cc(ie_cc,3)=sn_c(ie_cc,3)*r*r
      enddo
c
c***********************************
c
c finest mesh
c
      do ib=1,nbface
        read(5,*)iib,(kbface(ib,in),in=1,2)
      enddo
      do ib=1,nbface
        read(5,*)iib, (sb(ib,i),i=1,3)
      sb(ib,1)=sb(ib,1)*r
      sb(ib,2)=sb(ib,2)*r
      sb(ib,3)=sb(ib,3)*r*r
      enddo
      if(isfera.eq.1)then
      do ib = 1, nboun
        read(5,*)ii, ibpoin(ib,1), ibpoin(ib,2),irecog1(ib)
      enddo
      else
      do ib = 1, nboun
        read(5,*)ii, ibpoin(ib,1), ibpoin(ib,2)
      enddo
      endif
      do ib = 1, nboun
        read(5,*)ii, (sbb(ib,i),i=1,3)
      sbb(ib,1)=sbb(ib,1)*r
      sbb(ib,2)=sbb(ib,2)*r
      sbb(ib,3)=sbb(ib,3)*r*r
      enddo
      do ib = 1, nboun
        read(5,*)ii, (anp(ib,i),i=1,3)
      enddo
c
c coarse 1 mesh
c
      do ib_c=1,nbface_c
        read(35,*)iib_c,(kbface_c(ib_c,in),in=1,2)
      enddo
      do ib_c=1,nbface_c
        read(35,*)iib_c, (sb_c(ib_c,i),i=1,3)
      sb_c(ib_c,1)=sb_c(ib_c,1)*r
      sb_c(ib_c,2)=sb_c(ib_c,2)*r
      sb_c(ib_c,3)=sb_c(ib_c,3)*r*r
      enddo
      if(isfera.eq.1)then
      do ib_c = 1, nboun_c
        read(35,*)ii_c,ibpoin_c(ib_c,1),ibpoin_c(ib_c,2),irecog1_c(ib_c)
      enddo
      else
      do ib_c = 1, nboun_c
        read(35,*)ii_c, ibpoin_c(ib_c,1),ibpoin_c(ib_c,2)
      enddo
      endif
      do ib_c = 1, nboun_c
        read(35,*)ii_c, (sbb_c(ib_c,i),i=1,3)
      sbb_c(ib_c,1)=sbb_c(ib_c,1)*r
      sbb_c(ib_c,2)=sbb_c(ib_c,2)*r
      sbb_c(ib_c,3)=sbb_c(ib_c,3)*r*r
      enddo
      do ib_c = 1, nboun_c
        read(35,*)ii_c, (anp_c(ib_c,i),i=1,3)
      enddo
c
c coarse 2 mesh
c
      do ib_cc=1,nbface_cc
        read(45,*)iib_cc,(kbface_cc(ib_cc,in),in=1,2)
      enddo
      do ib_cc=1,nbface_cc
        read(45,*)iib_cc, (sb_cc(ib_cc,i),i=1,3)
      sb_cc(ib_cc,1)=sb_cc(ib_cc,1)*r
      sb_cc(ib_cc,2)=sb_cc(ib_cc,2)*r
      sb_cc(ib_cc,3)=sb_cc(ib_cc,3)*r*r
      enddo
      if(isfera.eq.1)then
      do ib_cc = 1, nboun_cc
        read(45,*)ii_cc,ibpoin_c(ib_cc,1),ibpoin_cc(ib_cc,2),
     + irecog1_cc(ib_cc)
      enddo
      else
      do ib_cc = 1, nboun_cc
        read(45,*)ii_cc, ibpoin_cc(ib_cc,1),ibpoin_cc(ib_cc,2)
      enddo
      endif
      do ib_cc = 1, nboun_cc
        read(45,*)ii_cc, (sbb_cc(ib_cc,i),i=1,3)
      sbb_cc(ib_cc,1)=sbb_cc(ib_cc,1)*r
      sbb_cc(ib_cc,2)=sbb_cc(ib_cc,2)*r
      sbb_cc(ib_cc,3)=sbb_cc(ib_cc,3)*r*r
      enddo
      do ib_cc = 1, nboun_cc
        read(45,*)ii_cc, (anp_cc(ib_cc,i),i=1,3)
      enddo

c
c finest mesh
c
500   format(1x,i10,3(1x,e26.20),2i10,3(1x,e26.20))
501   format(1x,i10,1x,e26.20,2i10,3(1x,e26.20))
      if(iper.eq.2)then
      read(5,*)ncoin
      do iic=1,ncoin
      read(5,*)icoin(iic,1),icoin(iic,2)
      enddo
      read(5,*)ncoin2
      do iic=1,ncoin2
      read(5,*)icoin2(iic,1),icoin2(iic,2)
      enddo
      read(5,*)ncoin4
      do iic=1,ncoin4
      read(5,*)icoin4(iic,1),icoin4(iic,2),icoin4(iic,3),icoin4(iic,4)
      enddo
      endif
      if(isfera.eq.1)then
      read(5,*)ncoin
      do iic=1,ncoin
      read(5,*)icoin(iic,1),icoin(iic,2)
      enddo
      endif
      !
      do ip=1,npoin
      timestp(ip)=delt
      enddo
      as=10.
c
c coarse mesh
c
      if(iper.eq.2)then
      read(35,*)ncoin_c
      do iic_c=1,ncoin_c
      read(35,*)icoin_c(iic_c,1),icoin_c(iic_c,2)
      enddo
      read(35,*)ncoin2_c
      do iic_c=1,ncoin2_c
      read(35,*)icoin2_c(iic_c,1),icoin2_c(iic_c,2)
      enddo
      read(35,*)ncoin4_c
      do iic_c=1,ncoin4_c
      read(35,*)icoin4_c(iic_c,1),icoin4_c(iic_c,2),
     + icoin4_c(iic_c,3),icoin4_c(iic_c,4)
      enddo
      endif
      if(isfera.eq.1)then
      read(35,*)ncoin_c
      do iic_c=1,ncoin_c
      read(35,*)icoin_c(iic_c,1),icoin(iic_c,2)
      enddo
      endif
c
c coarse 2 mesh
c
      if(iper.eq.2)then
      read(45,*)ncoin_cc
      do iic_cc=1,ncoin_cc
      read(45,*)icoin_cc(iic_cc,1),icoin_cc(iic_cc,2)
      enddo
      read(45,*)ncoin2_cc
      do iic_cc=1,ncoin2_cc
      read(45,*)icoin2_cc(iic_cc,1),icoin2_cc(iic_cc,2)
      enddo
      read(45,*)ncoin4_cc
      do iic_cc=1,ncoin4_cc
      read(45,*)icoin4_cc(iic_cc,1),icoin4_c(iic_cc,2),
     + icoin4_cc(iic_cc,3),icoin4_cc(iic_cc,4)
      enddo
      endif
      if(isfera.eq.1)then
      read(45,*)ncoin_cc
      do iic_cc=1,ncoin_cc
      read(45,*)icoin_cc(iic_cc,1),icoin_cc(iic_cc,2)
      enddo
      endif
c
c     call subroutines for the neighbours array
      write(6,*) "Neighbours definition"
      call neigh_definition(iedge,nedge,kbface,nbface,
     + count_neigh,neigh_list)
      call neigh_definition(iedge_c,nedge_c,kbface_c,nbface_c,
     + count_neigh_c,neigh_list_c)
      call neigh_definition(iedge_cc,nedge_cc,kbface_cc,nbface_cc,
     + count_neigh_cc,neigh_list_cc)
c
      !
      if (match_flag.eq.0) then
        write(6,*) "Matchings definition 1"
        call matching(coord,npoin,coord_c,npoin_c,match1)
        write(6,*) "Matching definition 2"
        call matching(coord_c,npoin_c,coord_cc,npoin_cc,match2)
        ! print in trial and trial2 files
*        do ip=1,npoin_c
*          write(7,*) match1(ip)
*        end do
*        do ip=1,npoin_cc
*          write(9,*) match2(ip)
*        end do
      else if (match_flag.eq.1) then
        do ip=1,npoin_c
          read(55,*) match1(ip)
        end do
        do ip=1,npoin_cc
          read(56,*) match2(ip)
        end do
      end if
c
c

 1000 format('m',2e17.5)
 2000 format('d',2e17.5)
      return
      end
c******************************************
c
c
c
      subroutine matching(pos_p,np_p,pos_c,np_c,match)!,unmatch)
c
c     The aim of this subroutine is to define a map between
c     points of coarse mesh and fine mesh. Match variable
c     is an array containing points of fine mesh that match with
c     points of the coarse in that particular position
c
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      dimension pos_p(npoinx,3), pos_c(npoinx,3)
      integer np_p,np_c,unmatch_count
      integer, dimension(npoinx) :: match!,unmatch
c
      tol=1.
      match=0
      !unmatch_count=0
      do ip=1,np_p
        ! define the position and check for matchings
        ! go out when the point is found
        x_p=pos_p(ip,1)
        y_p=pos_p(ip,2)
        z_p=pos_p(ip,3)
        do ic=1,np_c
          x_c=pos_c(ic,1)
          y_c=pos_c(ic,2)
          z_c=pos_c(ic,3)
          if ((abs(x_p-x_c).lt.tol).and.(abs(y_p-y_c).lt.tol)
     + .and.(abs(z_p-z_c).lt.tol)) then
              match(ic)=ip
              go to 99
          end if
        end do
       ! unmatch_count=unmatch_count+1
       ! unmatch(unmatch_count)=ip
  99    continue
      end do

      end subroutine
c
c
c
cc ----------------------------------------------------------
c
      subroutine neigh_definition(int_edg,nint_edg,ext_edg,next_edg,
     + count_nei,nei_list)
c
c     The aim of this subroutine is to compute the number of neighbours
c     for each point of the domain, and to individuate them.
c
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      dimension int_edg(nedgex,2)
      integer, dimension(nfacesx,1:2) :: ext_edg
      integer nint_edg,next_edg
      integer, dimension(npoinx) :: count_nei
      dimension nei_list(npoinx,10)
c
c
      count_nei = 0.
      nei_list = 0.
c
c     loop on internal edges
      do ie=1,nint_edg
         i1=int_edg(ie,1)
         i2=int_edg(ie,2)
         ! counting the neighbours
         count_nei(i1)=count_nei(i1)+1
         count_nei(i2)=count_nei(i2)+1
         ! storing them
         nei_list(i1,count_nei(i1))=i2
         nei_list(i2,count_nei(i2))=i1
      end do
c
c     loop on external edges
      do ie=1,next_edg
         i1=ext_edg(ie,1)
         i2=ext_edg(ie,2)
         ! counting the neighbours
         count_nei(i1)=count_nei(i1)+1
         count_nei(i2)=count_nei(i2)+1
         ! storing them
         nei_list(i1,count_nei(i1))=i2
         nei_list(i2,count_nei(i2))=i1
      end do
c
      end subroutine
c
c
c
      subroutine output(unk)
c
c***********************************************************
c writes final flowfield into results.d file
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      dimension unk(npoinx,4)
cc    change small domain
c     do ip = 1, npoin
c     alref2=128.0/15000.*100.
c       coord(ip,1) = coord(ip,1)*alref2
c       coord(ip,2) = coord(ip,2)*alref2
c       coord(ip,3) = coord(ip,3)*alref2
c     enddo !ip
c     writes values of the unknown
      write(11,*)1.
      write(11,*)npoin
c      write(11,'(a)')'point X Y u v rho press energy Mach'
      write(11,'(a)')'point X Y Z u v w temp Velocity'
      do ip=1,npoin
      u=unk(ip,1)
      v=unk(ip,2)
      w=unk(ip,3)
      vv=sqrt(u*u+v*v+w*w)
      th=unk(ip,4)
      scal=1.e+03
      write(11,100)ip,(coord(ip,i)/scal,i=1,3),u,v
     &,w,th,vv
      enddo
100   format(1x,i7,8e12.4)
      return
      end
c-----------
        subroutine meshstats
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      common/volume/vol(npoinx)

      vav=0.
      vmx=-1.e10
      vmn= 1.e10
      rav=0.
      rmx=-1.e10
      rmn= 1.e10
      do ip=1,npoin
      vav=vav+vol(ip)
      vmx=amax1(vmx,vol(ip))
      vmn=amin1(vmn,vol(ip))
      rav=rav+vol(ip)**(1./3.)
      rmx=amax1(rmx,vol(ip)**(1./3.))
      rmn=amin1(rmn,vol(ip)**(1./3.))
      enddo
      vav=vav/npoin
      rav=rav/npoin
      vsd=0.
      rsd=0.
      do ip=1,npoin
      vsd=vsd+(vol(ip)-vav)**2
      rsd=rsd+(vol(ip)**(1./3.)-rav)**2
      enddo
      vsd=sqrt(vsd/npoin)
      rsd=sqrt(rsd/npoin)
      write(6,205)vmx,vmn,vav,vsd,
     +            rmx,rmn,rav,rsd
  205 format(1x,'volmx,volmn,volav,volsd:',4e11.4/
     &       1x,'delmx,delmn,delav,delsd:',4e11.4)
        return
        end
        subroutine stats(unk)
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      dimension unk(npoinx,4)

      uav=0.
      umx=-1.e10
      umn= 1.e10
      vav=0.
      vmx=-1.e10
      vmn= 1.e10
      wav=0.
      wmx=-1.e10
      wmn= 1.e10
      tav=0.
      tmx=-1.e10
      tmn= 1.e10
      rav=0.
      do ip=1,npoin
      rav=rav+rh0(ip)
      uav=uav+rh0(ip)*unk(ip,1)
      umx=amax1(umx,unk(ip,1))
      umn=amin1(umn,unk(ip,1))
      vav=vav+rh0(ip)*unk(ip,2)
      vmx=amax1(vmx,unk(ip,2))
      vmn=amin1(vmn,unk(ip,2))
      wav=wav+rh0(ip)*unk(ip,3)
      wmx=amax1(wmx,unk(ip,3))
      wmn=amin1(wmn,unk(ip,3))
      tav=tav+rh0(ip)*unk(ip,4)
      tmx=amax1(tmx,unk(ip,4))
      tmn=amin1(tmn,unk(ip,4))
      enddo
      uav=uav/rav
      vav=vav/rav
      wav=wav/rav
      tav=tav/rav
      usd=0.
      vsd=0.
      wsd=0.
      tsd=0.
      do ip=1,npoin
      usd=usd+rh0(ip)*(unk(ip,1)-uav)**2
      vsd=vsd+rh0(ip)*(unk(ip,2)-vav)**2
      wsd=wsd+rh0(ip)*(unk(ip,3)-wav)**2
      tsd=tsd+rh0(ip)*(unk(ip,4)-tav)**2
      enddo
      usd=sqrt(usd/rav)
      vsd=sqrt(vsd/rav)
      wsd=sqrt(wsd/rav)
      tsd=sqrt(tsd/npoin)
      write(6,205)umx,umn,uav,usd,vmx,vmn,vav,vsd,
     +            wmx,wmn,wav,wsd,tmx,tmn,tav,tsd
  205 format(1x,'umx,umn,uav,usd:',4e11.4/
     &       1x,'vmx,vmn,vav,vsd:',4e11.4/
     &       1x,'wmx,wmn,wav,wsd:',4e11.4/
     &       1x,'tmx,tmn,tav,tsd:',4e11.4)

        return
        end

c---------------------------------------------
      subroutine absorbers(alpha)
c
c***********************************************************
c      Defines values of matrix alpha elements, ensuring smooth
c      decay of values prescribed on the boundary into computational
c      domain over specified in the input distances
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
c
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/xyz/coord(npoinx,3)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/absspec/dxabl,dxabr,dyabd,dyabu,dzabd,dzabu,
     &towxl,towxr,towyd,towyu,towzd,towzu,irl
      dimension alpha(npoinx)
c
      r=radious
      eps=1.e-07
c ----------------abers
      if ( irl.eq.1 ) then
      abxmax=-1.e+07
      abxmin=+1.e+07
      abymax=-1.e+07
      abymin=+1.e+07
      abzmax=-1.e+07
      abzmin=+1.e+07
      do ip=1,npoin
      abxmax=amax1(coord(ip,1),abxmax)
      abxmin=amin1(coord(ip,1),abxmin)
      abymax=amax1(coord(ip,2),abymax)
      abymin=amin1(coord(ip,2),abymin)
      abzmax=amax1(coord(ip,3),abzmax)
      abzmin=amin1(coord(ip,3),abzmin)
      enddo
      xabr=abxmax*r-dxabr
      xabl=abxmin*r+dxabl
      yabu=abymax*r-dyabu
      yabd=abymin*r+dyabd
      zabu=abzmax-dzabu
      zabd=abzmin+dzabd
c
      if(towxl.ne.0.)then
      toliL=1./towxL
      else
      tolil=0.
      endif
      if(towxr.ne.0.)then
      tolir=1./towxr
      else
      tolir=0.
      endif
      if(towyu.ne.0.)then
      toliu=1./towyu
      else
      toliu=0.
      endif
      if(towyd.ne.0.)then
      tolid=1./towyd
      else
      tolid=0.
      endif
      if(towzu.ne.0.)then
      toliuz=1./towzu
      else
      toliuz=0.
      endif
      if(towzd.ne.0.)then
      tolidz=1./towzd
      else
      tolidz=0.
      endif
       do ip=1,npoin
      relx=0.
      rely=0.
      relz=0.
       r0=coord(ip,1)*r
       r1=0.
       r2=0.
      relb=0.
      if(dxabl.ne.0.) r1=amax1(0., xabl-r0)/dxabL
      if(dxabr.ne.0.) r2=amax1(0., r0-xabr)/dxabR
      if((tolil+tolir).ne.0.)
     &relb=((toliL*r1)**2+(toliR*r2)**2)/(toliL*r1+toliR*r2+eps)
         relx=relb
       r1=0.
       r2=0.
      relb=0.
       r0=coord(ip,2)*r
         if(dyabd.ne.0.)r1=amax1(0.,-r0+yabd)/dyabd
         if(dyabu.ne.0.)r2=amax1(0., r0-yabu)/dyabu
      if((tolid+toliu).ne.0.)
     &relb=((tolid*r1)**2+(toliu*r2)**2)/(tolid*r1+toliu*r2+eps)
         rely=relb
       r1=0.
       r2=0.
      relb=0.
       r0=coord(ip,3)
         if(dzabd.ne.0.)r1=amax1(0.,-r0+zabd)/dzabd
         if(dzabu.ne.0.)r2=amax1(0., r0-zabu)/dzabu
      if((tolidz+toliuz).ne.0.)
     &relb=((tolidz*r1)**2+(toliuz*r2)**2)/(tolidz*r1+toliuz*r2+eps)
         relz=relb
c
      alpha1=(rely**2+relx**2)/(rely+relx+eps)
      alpha(ip)=(alpha1**2+relz**2)/(alpha1+relz+eps)
       enddo
      endif
c--------------------------------------end absorbers----------
      return
      end
c---------------------------------------------
c---------------------------------
      subroutine mpdatam(uuu,unk,iflag)
cc      implicit double precision(a-h,o-z)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      dimension afreelc(5)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz/coord(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume/vol(npoinx)
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/terms/sbb(nfacesx,3)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension aun(nedgex),aunb(nfacesx),iiflg(npoinx)
      dimension aunbb1(nfacesx)
      dimension uuu(npoinx,3)
      dimension divv(npoinx)
      dimension dpdxa(npoinx),dpdya(npoinx),dpdza(npoinx)
      dimension unkmx(npoinx),unkmn(npoinx),rhin(npoinx),rhout(npoinx)
      dimension cp(npoinx),cn(npoinx),switch(npoinx),switch1(npoinx)
      dimension arh2(npoinx),unk(npoinx)
      partp(yy)= amax1(0.,yy)
      partn(yy)=-amin1(0.,yy)
c
c
      iord=2
      npass=iord
      do ip = 1,npoin
      divv(ip)=0.
      switch1(ip)=switch(ip)
      switch(ip)=1.
      switch1(ip)=1.
      iiflg(ip)=0.
      enddo
      do ied = 1,nedge
      aun(ied)=0.
      enddo
      do ied = 1,nbface
      aunb(ied)=0.
      enddo
      do ied=1,nboun
      aunbb1(ied)=0.
      enddo
c
      ifct=1
c negative number of passes
      npass=abs(npass)
      tol=1.e-10
c
c loop over mpdata passes
      do ip=1,npoin
      unkmx(ip)=-1.e15
      unkmn(ip)= 1.e15
      enddo
                 do 3 ipass=1,npass+1
c
c
      do ip = 1,npoin
      arh2(ip)=0.
      if(ipass.eq.1)switch(ip)=1.
      if(ipass.ne.1)switch(ip)=switch1(ip)
      enddo
c  -----internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      if(ipass.eq.1)then
      aun(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
      fl=unk(i1)*amax1(0.,aun(ied))+unk(i2)*amin1(0.,aun(ied))
      else
      fl=aun(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aun(ied)
      divv(i2)=divv(i2)-aun(ied)
      enddo
c  ----- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      if(ipass.eq.1)then
      aunb(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
c
      fl=unk(i1)*amax1(0.,aunb(ied))+unk(i2)*amin1(0.,aunb(ied))
      else
      fl=aunb(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aunb(ied)
      divv(i2)=divv(i2)-aunb(ied)
      enddo
c  ----- boundary points
      do  ied = 1,nboun
      i1=ibpoin(ied,1)
      sx=sbb(ied,1)
      sy=sbb(ied,2)
      sz=sbb(ied,3)
      x=coord(i1,1)
      y=coord(i1,2)
      z=coord(i1,3)
c
      if(ipass.eq.1)then
      aunbb1(ied)=-divv(i1)
      afreelc(iflag)=afree(iflag)
      afreelc(1)=ue(i1)
      fl1=afreelc(iflag)*amin1(0.,aunbb1(ied))
     &   +unk(i1)*amax1(0.,aunbb1(ied))
      else
      fl1=0.
      endif
      if(ibpoin(ied,2).eq.11
     &.or.ibpoin(ied,2).eq.5
     &.or.ibpoin(ied,2).eq.6
     &)then
      fl1=0.
      aunbb1(ied)=0.
      endif
      arh2(i1)=arh2(i1)+fl1*switch(i1)
      enddo
c---------
c update unknowns
      do ip = 1,npoin
      unk(ip)=unk(ip)-arh2(ip)*timestp(ip)/vol(ip)/rh0(ip)
      enddo
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
c---------
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(ipass.eq.iord) goto 6
c
c NORMAL VELOCITY FOR THE SECOND PASS
c
      call derivative(unk,dpdxa,1)
      call derivative(unk,dpdya,2)
      call derivative(unk,dpdza,3)
c
c ----  internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      un=abs(aun(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aun(ied)= un-aun(ied)*term*timestp(1)*0.5
      enddo
c ---- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      un=abs(aunb(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aunb(ied)= un-aunb(ied)*term*timestp(1)*0.5
      enddo
      do  ied = 1,nboun
      aunbb1(ied)=0.
      enddo
c  non-isscilatory option
                         if(ifct.eq.1)then
cccc      goto 1000
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
      enddo
      do  ib = 1,nbface
      i1=kbface(ib,1)
      i2=kbface(ib,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
        enddo
cccc1000  continue

      do ip = 1,npoin
      rhin(ip)=0.
      rhout(ip)=0.
      cp(ip)=0.
      cn(ip)=0.
      enddo
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      rhin(i1)=rhin(i1)+partn(aun(ied))
      rhout(i1)=rhout(i1)+partp(aun(ied))
      rhin(i2)=rhin(i2)+partp(aun(ied))
      rhout(i2)=rhout(i2)+partn(aun(ied))
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      rhin(i1)=rhin(i1)+partn(aunb(ied))
      rhout(i1)=rhout(i1)+partp(aunb(ied))
      rhin(i2)=rhin(i2)+partp(aunb(ied))
      rhout(i2)=rhout(i2)+partn(aunb(ied))
      enddo
c
      do ip=1,npoin
      cp(ip)=(unkmx(ip)-unk(ip))*vol(ip)*rh0(ip)/
     &(rhin(ip)*timestp(ip)+tol)
      cn(ip)=(unk(ip)-unkmn(ip))*vol(ip)*rh0(ip)/
     &(rhout(ip)*timestp(ip)+tol)
      enddo
c
c limited antidiffusive  velocities:
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      if(aun(ied).gt.0.)then
      aun(ied)=aun(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aun(ied)=aun(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      if(aunb(ied).gt.0.)then
      aunb(ied)=aunb(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aunb(ied)=aunb(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
                                    endif
c end of non-isscilatory option
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c end of loop over mpdata passes
  3                continue
  6                continue

      return
      end
c******************************************
c---------------------------------
      subroutine mpdatams(uuu,unk,iflag)
cc      implicit double precision(a-h,o-z)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz/coord(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume/vol(npoinx)
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/terms/sbb(nfacesx,3)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension aun(nedgex),aunb(nfacesx),iiflg(npoinx)
      dimension uuu(npoinx,3)
      dimension divv(npoinx)
      dimension dpdxa(npoinx),dpdya(npoinx),dpdza(npoinx)
      dimension unkmx(npoinx),unkmn(npoinx),rhin(npoinx),rhout(npoinx)
      dimension cp(npoinx),cn(npoinx),switch(npoinx),switch1(npoinx)
      dimension arh2(npoinx),unk(npoinx),aunbb1(nfacesx)
      partp(yy)= amax1(0.,yy)
      partn(yy)=-amin1(0.,yy)
c
c
      iord=2
      npass=iord
      do ip = 1,npoin
      divv(ip)=0.
      switch1(ip)=switch(ip)
      switch(ip)=1.
      switch1(ip)=1.
      iiflg(ip)=0.
      enddo
      do ied = 1,nedge
      aun(ied)=0.
      enddo
      do ied = 1,nbface
      aunb(ied)=0.
      enddo
      do ied = 1,nboun
      aunbb1(ied)=0.
      enddo
c
      ifct=1
c negative number of passes
      npass=abs(npass)
      tol=1.e-10
c
c loop over mpdata passes
      do ip=1,npoin
      unkmx(ip)=-1.e15
      unkmn(ip)= 1.e15
      enddo
                 do 3 ipass=1,npass+1
c
c
      do ip = 1,npoin
      arh2(ip)=0.
      if(ipass.eq.1)switch(ip)=1.
      if(ipass.ne.1)switch(ip)=switch1(ip)
      enddo
c  -----internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      if(ipass.eq.1)then
      aun(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
      fl=unk(i1)*amax1(0.,aun(ied))+unk(i2)*amin1(0.,aun(ied))
      else
      fl=aun(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aun(ied)
      divv(i2)=divv(i2)-aun(ied)
      enddo
c  ----- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      if(ipass.eq.1)then
      aunb(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
c
      fl=unk(i1)*amax1(0.,aunb(ied))+unk(i2)*amin1(0.,aunb(ied))
      else
      fl=aunb(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aunb(ied)
      divv(i2)=divv(i2)-aunb(ied)
      enddo
c  ----- boundary points
      do  ied = 1,nboun
      i1=ibpoin(ied,1)
c
      aunbb1(ied)=-divv(i1)
      fl1=0.
      arh2(i1)=arh2(i1)+fl1*switch(i1)
      enddo
c---------
c---------
      call cyclicr(arh2)
c update unknowns
      do ip = 1,npoin
      unk(ip)=unk(ip)-arh2(ip)*timestp(ip)/vol(ip)/rh0(ip)
      enddo
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
c---------
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(ipass.eq.iord) goto 6
c
c NORMAL VELOCITY FOR THE SECOND PASS
c
      if(iflag.lt.3)then
      call derivative(unk,dpdxa,1)
      call derivativeminus(unk,dpdya,2)
      call derivative(unk,dpdza,3)
      else
      call derivative(unk,dpdxa,1)
      call derivative(unk,dpdya,2)
      call derivative(unk,dpdza,3)
      endif
c
c ----  internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      un=abs(aun(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aun(ied)= un-aun(ied)*term*timestp(1)*0.5
      enddo
c ---- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      un=abs(aunb(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aunb(ied)= un-aunb(ied)*term*timestp(1)*0.5
      enddo
c  non-isscilatory option
                         if(ifct.eq.1)then
cccc      goto 1000
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
      enddo
      do  ib = 1,nbface
      i1=kbface(ib,1)
      i2=kbface(ib,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
        enddo

ccccccccccccccccccccccccccccccc

      do ib = 1,nboun
      if(ibpoin(ib,2).eq.5)then
      i1=ibpoin(ib,1)
      p1=unk(i1)
      ir11=irecog1(ib)
      pp1=unk(ir11)
      if(iflag.eq.1.or.iflag.eq.2)then
      pp1=-pp1
      endif
      unkmx(i1)=amax1(unkmx(i1),p1,pp1)
      unkmn(i1)=amin1(unkmn(i1),p1,pp1)
      endif
        enddo
      do ic=1,ncoin
      i1=icoin(ic,1)
      i2=icoin(ic,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2,unkmx(i2))
      unkmn(i1)=amin1(unkmn(i1),p1,p2,unkmn(i2))
      unkmx(i2)=amax1(unkmx(i2),p1,p2,unkmx(i1))
      unkmn(i2)=amin1(unkmn(i2),p1,p2,unkmn(i1))
      enddo


ccccccccccccccccccccccccccccccccccccccc

cccc1000  continue

      do ip = 1,npoin
      rhin(ip)=0.
      rhout(ip)=0.
      cp(ip)=0.
      cn(ip)=0.
      enddo
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      rhin(i1)=rhin(i1)+partn(aun(ied))
      rhout(i1)=rhout(i1)+partp(aun(ied))
      rhin(i2)=rhin(i2)+partp(aun(ied))
      rhout(i2)=rhout(i2)+partn(aun(ied))
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      rhin(i1)=rhin(i1)+partn(aunb(ied))
      rhout(i1)=rhout(i1)+partp(aunb(ied))
      rhin(i2)=rhin(i2)+partp(aunb(ied))
      rhout(i2)=rhout(i2)+partn(aunb(ied))
      enddo
c boundary edges closure
      call cyclicr(rhin)
      call cyclicr(rhout)
c
      do ip=1,npoin
      cp(ip)=(unkmx(ip)-unk(ip))*vol(ip)*rh0(ip)/
     &(rhin(ip)*timestp(ip)+tol)
      cn(ip)=(unk(ip)-unkmn(ip))*vol(ip)*rh0(ip)/
     &(rhout(ip)*timestp(ip)+tol)
      enddo
c
c limited antidiffusive  velocities:
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      if(aun(ied).gt.0.)then
      aun(ied)=aun(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aun(ied)=aun(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      if(aunb(ied).gt.0.)then
      aunb(ied)=aunb(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aunb(ied)=aunb(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
                                    endif
c end of non-isscilatory option
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c end of loop over mpdata passes
  3                continue
  6                continue


      return
      end
c******************************************
c---------------------------------
      subroutine mpdatamper(uuu,unk,iflag)
cc      implicit double precision(a-h,o-z)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
       parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per2/  icoin2(nfacesx,2),ncoin2
      common/per4/  icoin4(nfacesx,4),ncoin4
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/xyz/coord(npoinx,3)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/volume/vol(npoinx)
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/terms/sbb(nfacesx,3)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension aun(nedgex),aunb(nfacesx),iiflg(npoinx)
      dimension uuu(npoinx,3)
      dimension divv(npoinx)
      dimension dpdxa(npoinx),dpdya(npoinx),dpdza(npoinx)
      dimension unkmx(npoinx),unkmn(npoinx),rhin(npoinx),rhout(npoinx)
      dimension cp(npoinx),cn(npoinx),switch(npoinx),switch1(npoinx)
      dimension arh2(npoinx),unk(npoinx),aunbb1(nfacesx)
      partp(yy)= amax1(0.,yy)
      partn(yy)=-amin1(0.,yy)
c
c
      iord=2
      npass=iord
      do ip = 1,npoin
      divv(ip)=0.
      switch1(ip)=switch(ip)
      switch(ip)=1.
      switch1(ip)=1.
      iiflg(ip)=0.
      enddo
      do ied = 1,nedge
      aun(ied)=0.
      enddo
      do ied = 1,nbface
      aunb(ied)=0.
      enddo
      do ied = 1,nboun
      aunbb1(ied)=0.
      enddo
c
      ifct=1
c negative number of passes
      npass=abs(npass)
      tol=1.e-10
c
c loop over mpdata passes
      do ip=1,npoin
      unkmx(ip)=-1.e15
      unkmn(ip)= 1.e15
      enddo
                 do 3 ipass=1,npass+1
c
c
      do ip = 1,npoin
      arh2(ip)=0.
      if(ipass.eq.1)switch(ip)=1.
      if(ipass.ne.1)switch(ip)=switch1(ip)
      enddo
c  -----internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      if(ipass.eq.1)then
      aun(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
      fl=unk(i1)*amax1(0.,aun(ied))+unk(i2)*amin1(0.,aun(ied))
      else
      fl=aun(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aun(ied)
      divv(i2)=divv(i2)-aun(ied)
      enddo
c  ----- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      if(ipass.eq.1)then
      aunb(ied)=
     & (uuu(i1,1)+uuu(i2,1))*0.5*sx
     &+(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &+(uuu(i1,3)+uuu(i2,3))*0.5*sz
c
      fl=unk(i1)*amax1(0.,aunb(ied))+unk(i2)*amin1(0.,aunb(ied))
      else
      fl=aunb(ied)
      endif
      arh2(i1)=arh2(i1)+fl*(switch(i1)+switch(i2))*0.5
      arh2(i2)=arh2(i2)-fl*(switch(i1)+switch(i2))*0.5
      divv(i1)=divv(i1)+aunb(ied)
      divv(i2)=divv(i2)-aunb(ied)
      enddo
c  ----- boundary points
      do  ied = 1,nboun
      i1=ibpoin(ied,1)
c
      aunbb1(ied)=0.
      fl1=0.
      arh2(i1)=arh2(i1)+fl1*switch(i1)
      enddo
c---------
c---------
      call cyclicr(arh2)
        call cyclicr2(arh2)
        call cyclicr4(arh2)
c update unknowns
      do ip = 1,npoin
      unk(ip)=unk(ip)-arh2(ip)*timestp(ip)/vol(ip)/rh0(ip)
      enddo
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
c---------
c ~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(ipass.eq.iord) goto 6
c
c NORMAL VELOCITY FOR THE SECOND PASS
c
      call derivative(unk,dpdxa,1)
      call derivative(unk,dpdya,2)
      call derivative(unk,dpdza,3)
c
c ----  internal edges
      do  ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
c
      un=abs(aun(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aun(ied)= un-aun(ied)*term*timestp(1)*0.5
      enddo
c ---- semi edges
      do  ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
c
      un=abs(aunb(ied))*(unk(i2)-unk(i1))*0.5
       volt=0.5*(vol(i1)*rh0(i1)+vol(i2)*rh0(i2))
        dpdx=0.5*(dpdxa(i1)+dpdxa(i2))/volt
        dpdy=0.5*(dpdya(i1)+dpdya(i2))/volt
        dpdz=0.5*(dpdza(i1)+dpdza(i2))/volt
      term=   (uuu(i1,1)+uuu(i2,1))*0.5*dpdx
     &       +(uuu(i1,2)+uuu(i2,2))*0.5*dpdy
     &       +(uuu(i1,3)+uuu(i2,3))*0.5*dpdz
      aunb(ied)= un-aunb(ied)*term*timestp(1)*0.5
      enddo
c  non-isscilatory option
                         if(ifct.eq.1)then
cccc      goto 1000
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
      enddo
      do  ib = 1,nbface
      i1=kbface(ib,1)
      i2=kbface(ib,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2)
      unkmn(i1)=amin1(unkmn(i1),p1,p2)
      unkmx(i2)=amax1(unkmx(i2),p2,p1)
      unkmn(i2)=amin1(unkmn(i2),p2,p1)
        enddo

ccccccccccccccccccccccccccccccc

      do ic=1,ncoin
      i1=icoin(ic,1)
      i2=icoin(ic,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2,unkmx(i2))
      unkmn(i1)=amin1(unkmn(i1),p1,p2,unkmn(i2))
      unkmx(i2)=amax1(unkmx(i2),p1,p2,unkmx(i1))
      unkmn(i2)=amin1(unkmn(i2),p1,p2,unkmn(i1))
      enddo
      do ic=1,ncoin2
      i1=icoin2(ic,1)
      i2=icoin2(ic,2)
      p1=unk(i1)
      p2=unk(i2)
      unkmx(i1)=amax1(unkmx(i1),p1,p2,unkmx(i2))
      unkmn(i1)=amin1(unkmn(i1),p1,p2,unkmn(i2))
      unkmx(i2)=amax1(unkmx(i2),p1,p2,unkmx(i1))
      unkmn(i2)=amin1(unkmn(i2),p1,p2,unkmn(i1))
      enddo
      do ic=1,ncoin4
      i1=icoin4(ic,1)
      i2=icoin4(ic,2)
      i3=icoin4(ic,3)
      i4=icoin4(ic,4)
      p1=unk(i1)
      p1=unk(i1)
      p2=unk(i2)
      p3=unk(i3)
      p4=unk(i4)
      unkmx(i1)=amax1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmn(i1)=amin1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmx(i2)=amax1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmn(i2)=amin1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmx(i3)=amax1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmn(i3)=amin1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmx(i4)=amax1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      unkmn(i4)=amin1(unkmx(i1),unkmx(i2),unkmx(i3),unkmx(i4)
     &,p1,p2,p3,p4)
      enddo
ccccccccccccccccccccccccccccccccccccccc

cccc1000  continue

      do ip = 1,npoin
      rhin(ip)=0.
      rhout(ip)=0.
      cp(ip)=0.
      cn(ip)=0.
      enddo
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      rhin(i1)=rhin(i1)+partn(aun(ied))
      rhout(i1)=rhout(i1)+partp(aun(ied))
      rhin(i2)=rhin(i2)+partp(aun(ied))
      rhout(i2)=rhout(i2)+partn(aun(ied))
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      rhin(i1)=rhin(i1)+partn(aunb(ied))
      rhout(i1)=rhout(i1)+partp(aunb(ied))
      rhin(i2)=rhin(i2)+partp(aunb(ied))
      rhout(i2)=rhout(i2)+partn(aunb(ied))
      enddo
c boundary edges closure
      call cyclicr(rhin)
      call cyclicr(rhout)
      call cyclicr2(rhin)
      call cyclicr2(rhout)
      call cyclicr4(rhin)
      call cyclicr4(rhout)
c
      do ip=1,npoin
      cp(ip)=(unkmx(ip)-unk(ip))*vol(ip)*rh0(ip)/
     &(rhin(ip)*timestp(ip)+tol)
      cn(ip)=(unk(ip)-unkmn(ip))*vol(ip)*rh0(ip)/
     &(rhout(ip)*timestp(ip)+tol)
      enddo
c
c limited antidiffusive  velocities:
c internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      if(aun(ied).gt.0.)then
      aun(ied)=aun(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aun(ied)=aun(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
c boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      if(aunb(ied).gt.0.)then
      aunb(ied)=aunb(ied)*amin1(1.,cp(i2),cn(i1))
      else
      aunb(ied)=aunb(ied)*amin1(1.,cn(i2),cp(i1))
      endif
      enddo
                                    endif
c end of non-isscilatory option
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c end of loop over mpdata passes
  3                continue
  6                continue


      return
      end
cc******************************************
      subroutine derivative(unk,dpdxa,idx)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension dpdxa(npoinx)
      dimension unk(npoinx)
c
c
      tol=1.e-07
c  derivatives
      do ip = 1,npoin
      dpdxa(ip)=0.
      enddo
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      if(i1.eq.1.or.i2.eq.1)then
      a=1.
      endif
      sx=sb(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo

      IF(ISFERA.EQ.1) THEN
      do ied=1,nboun
      i1=ibpoin(ied,1)
      if(ibpoin(ied,2).eq.5.and.idx.eq.2)then
      ir11=irecog1(ied)
      if(ir11.eq.0)then
      write(6,*)'unmatched point'
      endif
      sx=sbb(ied,idx)
      u=(unk(ir11)+unk(i1))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      if(ibpoin(ied,2).eq.11.and.idx.eq.3)then
      sx=sbb(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      enddo
      call cyclicr(dpdxa)
      Endif
      if(isfera.eq.0.and.iper.eq.0)then
        do ied=1,nboun
        i1=ibpoin(ied,1)
      sx=sbb(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      enddo
      ENDIF
      if(iper.eq.2)then
      call cyclicr(dpdxa)
      call cyclicr2(dpdxa)
      call cyclicr4(dpdxa)
      Endif

      return
      end
c-----------
      subroutine derivative_c(unk,dpdxa,idx)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/coin_c/ncoin_c,icoin_c(nfacesx*2,2),irecog1_c(nfacesx)
      common/bpoints_c/ibpoin_c(nfacesx,2),anp_c(nfacesx,3)
      common/volume_c/vol_c(npoinx)
      common/terms_c/sbb_c(nfacesx,3)
      common/edg_c/iedge_c(nedgex,2)
      common/n_c/npoin_c,nedge_c,nbface_c,nbpoin_c,nboun_c
      common/wge_c/sn_c(nedgex,3)
      common/bound_c/kbface_c(nfacesx,1:2),sb_c(nfacesx,3)
      dimension dpdxa(npoinx)
      dimension unk(npoinx)
c
c
      tol=1.e-07
c  derivatives
      do ip = 1,npoin_c
      dpdxa(ip)=0.
      enddo
      do ied = 1,nedge_c
      i1=iedge_c(ied,1)
      i2=iedge_c(ied,2)
      sx=sn_c(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied = 1,nbface_c
      i1=kbface_c(ied,1)
      i2=kbface_c(ied,2)
      if(i1.eq.1.or.i2.eq.1)then
      a=1.
      endif
      sx=sb_c(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo

      IF(ISFERA.EQ.1) THEN
      do ied=1,nboun_c
      i1=ibpoin_c(ied,1)
      if(ibpoin_c(ied,2).eq.5.and.idx.eq.2)then
      ir11=irecog1_c(ied)
      if(ir11.eq.0)then
      write(6,*)'unmatched point'
      endif
      sx=sbb_c(ied,idx)
      u=(unk(ir11)+unk(i1))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      if(ibpoin_c(ied,2).eq.11.and.idx.eq.3)then
      sx=sbb_c(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      enddo
      call cyclicr(dpdxa)
      Endif
      if(isfera.eq.0.and.iper.eq.0)then
        do ied=1,nboun_c
        i1=ibpoin_c(ied,1)
      sx=sbb_c(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      enddo
      ENDIF
      if(iper.eq.2)then
      call cyclicr(dpdxa)
      call cyclicr2(dpdxa)
      call cyclicr4(dpdxa)
      Endif

      return
      end
c
cc******************************************
      subroutine derivative_cc(unk,dpdxa,idx)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/per/iper
      common/sf/isfera
      common/coin_cc/ncoin_cc,icoin_cc(nfacesx*2,2),irecog1_cc(nfacesx)
      common/bpoints_cc/ibpoin_cc(nfacesx,2),anp_cc(nfacesx,3)
      common/volume_cc/vol_cc(npoinx)
      common/terms_cc/sbb_cc(nfacesx,3)
      common/edg_cc/iedge_cc(nedgex,2)
      common/n_cc/npoin_cc,nedge_cc,nbface_cc,nbpoin_cc,nboun_cc
      common/wge_cc/sn_cc(nedgex,3)
      common/bound_cc/kbface_cc(nfacesx,1:2),sb_cc(nfacesx,3)
      dimension dpdxa(npoinx)
      dimension unk(npoinx)
c
c
      tol=1.e-07
c  derivatives
      do ip = 1,npoin_cc
      dpdxa(ip)=0.
      enddo
      do ied = 1,nedge_cc
      i1=iedge_cc(ied,1)
      i2=iedge_cc(ied,2)
      sx=sn_cc(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied = 1,nbface_cc
      i1=kbface_cc(ied,1)
      i2=kbface_cc(ied,2)
      if(i1.eq.1.or.i2.eq.1)then
      a=1.
      endif
      sx=sb_cc(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo

      IF(ISFERA.EQ.1) THEN
      do ied=1,nboun_cc
      i1=ibpoin_cc(ied,1)
      if(ibpoin_cc(ied,2).eq.5.and.idx.eq.2)then
      ir11=irecog1_cc(ied)
      if(ir11.eq.0)then
      write(6,*)'unmatched point'
      endif
      sx=sbb_cc(ied,idx)
      u=(unk(ir11)+unk(i1))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      if(ibpoin_cc(ied,2).eq.11.and.idx.eq.3)then
      sx=sbb_cc(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      enddo
      call cyclicr(dpdxa)
      Endif
      if(isfera.eq.0.and.iper.eq.0)then
        do ied=1,nboun_cc
        i1=ibpoin_cc(ied,1)
      sx=sbb_cc(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      enddo
      ENDIF
      if(iper.eq.2)then
      call cyclicr(dpdxa)
      call cyclicr2(dpdxa)
      call cyclicr4(dpdxa)
      Endif

      return
      end
c-----------
      subroutine derivativeminus(unk,dpdxa,idx)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension dpdxa(npoinx)
      dimension unk(npoinx)
c
c
      tol=1.e-07
c  derivatives
      do ip = 1,npoin
      dpdxa(ip)=0.
      enddo
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied=1,nboun
      i1=ibpoin(ied,1)
      if(ibpoin(ied,2).eq.5.and.idx.eq.2)then
      ir11=irecog1(ied)
      sx=sbb(ied,idx)
      u=(-unk(ir11)+unk(i1))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      if(ibpoin(ied,2).eq.11.and.idx.eq.3)then
      sx=sbb(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      enddo
      call cyclicr(dpdxa)
      return
      end
c-----------
c-----------
      subroutine derivativediv(unk,dpdxa,idx)
c rhight hand side - contribution from inviscid fluxes
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension dpdxa(npoinx)
      dimension unk(npoinx)
c
c
      tol=1.e-07
c  derivatives
      do ip = 1,npoin
      dpdxa(ip)=0.
      enddo
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,idx)
      u=(unk(i1)+unk(i2))*0.5
      dpdxa(i1)=dpdxa(i1)+sx*u
      dpdxa(i2)=dpdxa(i2)-sx*u
      enddo
      do ied=1,nboun
      i1=ibpoin(ied,1)
      if(ibpoin(ied,2).eq.5.and.idx.eq.2)then
      u=0.
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      if(ibpoin(ied,2).eq.11.and.idx.eq.3)then
      sx=sbb(ied,idx)
      u=unk(i1)
      dpdxa(i1)=dpdxa(i1)+sx*u
      endif
      enddo
      call cyclicr(dpdxa)
      return
      end
c-----------
cc---------------------
      subroutine timecourant(uuu)
c***********************************************************
c provides addmissible by CFL condition timesteps   at every
c mesh point
c***********************************************************
cc      implicit double precision(a-h,o-z)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sss/delt,adapt
      common/cflpass/cfl
      common/tsteps/timestp(npoinx),istst,mrestart,time
      common/volume/vol(npoinx)
      common/terms/sbb(nfacesx,3)
      common/bpoints/ibpoin(nfacesx,2),anp(nfacesx,3)
      common/edg/iedge(nedgex,2)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/wge/sn(nedgex,3)
      common/bound/kbface(nfacesx,1:2),sb(nfacesx,3)
      dimension uuu(npoinx,3),cour(npoinx)
c
      do ip = 1,npoin
      cour(ip)=0.
      enddo
c  loop over the internal edges
      do ied = 1,nedge
      i1=iedge(ied,1)
      i2=iedge(ied,2)
      sx=sn(ied,1)
      sy=sn(ied,2)
      sz=sn(ied,3)
      aun=(uuu(i1,1)+uuu(i2,1))*0.5*sx
     &   +(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &   +(uuu(i1,3)+uuu(i2,3))*0.5*sz
      cour(i1)=cour(i1)+abs(aun)
      cour(i2)=cour(i2)+abs(aun)
      enddo
c
c  loop over boundary semi edges
      do ied = 1,nbface
      i1=kbface(ied,1)
      i2=kbface(ied,2)
      sx=sb(ied,1)
      sy=sb(ied,2)
      sz=sb(ied,3)
      aun=(uuu(i1,1)+uuu(i2,1))*0.5*sx
     &   +(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &   +(uuu(i1,3)+uuu(i2,3))*0.5*sz
      cour(i1)=cour(i1)+abs(aun)
      cour(i2)=cour(i2)+abs(aun)
      enddo
c  loop over boundary faces edges
      do ied=1,nboun
      i1=ibpoin(ied,1)
      sx=sbb(ied,1)
      sy=sbb(ied,2)
      sz=sbb(ied,3)
      aun=(uuu(i1,1)+uuu(i2,1))*0.5*sx
     &   +(uuu(i1,2)+uuu(i2,2))*0.5*sy
     &   +(uuu(i1,3)+uuu(i2,3))*0.5*sz
      cour(i1)=cour(i1)+abs(aun)
      enddo
      do ip = 1,npoin
      timestp(ip)=cfl*vol(ip)/cour(ip)
      enddo
      if(istst.ne.1)then
      tglob=1.e+15
      do ip=1,npoin
      tglob=amin1(tglob,timestp(ip))
      enddo
      if(istst.eq.2)tglob=delt
      do ip=1,npoin
      timestp(ip)=tglob
      enddo
      endif

      return
      end
cc---------------------
      subroutine cyclicr(unk)
       parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/xyz/coord(npoinx,3)
      common/poles/nsouth,north,isouth(nfacesx),inorth(nfacesx)
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/volume/vol(npoinx)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      dimension unk(npoinx),uc(nfacesx)
      do ic=1,ncoin
      il=icoin(ic,1)
      ir=icoin(ic,2)
      uc(ic)=unk(il)+unk(ir)
      enddo
      do ic=1,ncoin
      il=icoin(ic,1)
      ir=icoin(ic,2)
      unk(il)=uc(ic)
      unk(ir)=uc(ic)
      enddo

      return
      end
c------------------------
      subroutine cyclicr2(perd)

      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/per2/  icoin2(nfacesx,2),ncoin2
      dimension  perd(npoinx),uc(nfacesx)

      do ic = 1, ncoin2
        il = icoin2(ic,1)
        ir = icoin2(ic,2)
        uc(ic) = perd(il) + perd(ir)
      enddo
      do ic = 1, ncoin2
        il = icoin2(ic,1)
        ir = icoin2(ic,2)
        perd(il) = uc(ic)
        perd(ir) = uc(ic)
      enddo
      return
      end
c---------
c---------
      subroutine cyclicr4(perd)

      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/volume/vol(npoinx)
      common/per4/  icoin4(nfacesx,4),ncoin4
      dimension  perd(npoinx),uc(nfacesx)
      do ic = 1, ncoin4
        i1 = icoin4(ic,1)
        i2 = icoin4(ic,2)
        i3 = icoin4(ic,3)
        i4 = icoin4(ic,4)
      uc(ic)=perd(i1)+perd(i2)+perd(i3)+perd(i4)
      enddo
      do ic = 1, ncoin4
        i1 = icoin4(ic,1)
        i2 = icoin4(ic,2)
        i3 = icoin4(ic,3)
        i4 = icoin4(ic,4)
        perd(i1) = uc(ic)
        perd(i2) = uc(ic)
        perd(i3) = uc(ic)
        perd(i4) = uc(ic)
      enddo
      return
      end
c---------------------------------------------
      subroutine velprd(unk,frc,uuu,uuu1,unkdash)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
c
      common/sf/isfera
      common/metric/cosa(npoinx),gmm(npoinx),radious
      common/volume/vol(npoinx)
      common/xyz/coord(npoinx,3)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/tsteps/timestp(npoinx),istst,mrestart,time
      dimension unkdash(npoinx,3),ff1(npoinx),ff2(npoinx),ff3(npoinx)
      dimension uuu(npoinx,3),uuu1(npoinx,3)
      dimension unk(npoinx,4),frc(npoinx,4)
      dimension pprpot(npoinx)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps

      do ip=1,npoin
      unkdash(ip,1)=unk(ip,1)+timestp(ip)*frc(ip,1)
      unkdash(ip,2)=unk(ip,2)+timestp(ip)*frc(ip,2)
      unkdash(ip,3)=unk(ip,3)+timestp(ip)*frc(ip,3)
      enddo
      call derivative(unkdash(1,1),ff1,1)
      if(isfera.eq.1) call derivativeminus(unkdash(1,1),ff2,2)
      if(isfera.eq.0) call derivative(unkdash(1,1),ff2,2)
      call derivative(unkdash(1,1),ff3,3)
      do ip=1,npoin
      ff1(ip)=ff1(ip)/vol(ip)
      ff2(ip)=ff2(ip)/vol(ip)
      ff3(ip)=ff3(ip)/vol(ip)
      enddo
      do ip=1,npoin
      uuu1(ip,1)=unkdash(ip,1)-timestp(ip)*
     &(uuu(ip,1)*ff1(ip)+uuu(ip,2)*ff2(ip)+uuu(ip,3)*ff3(ip))
      enddo
      call derivative(unkdash(1,2),ff1,1)
      if(isfera.eq.1) call derivativeminus(unkdash(1,2),ff2,2)
      if(isfera.eq.0) call derivative(unkdash(1,2),ff2,2)
      call derivative(unkdash(1,2),ff3,3)
      do ip=1,npoin
      ff1(ip)=ff1(ip)/vol(ip)
      ff2(ip)=ff2(ip)/vol(ip)
      ff3(ip)=ff3(ip)/vol(ip)
      enddo
      do ip=1,npoin
      uuu1(ip,2)=unkdash(ip,2)-timestp(ip)*
     &(uuu(ip,1)*ff1(ip)+uuu(ip,2)*ff2(ip)+uuu(ip,3)*ff3(ip))
      enddo
      call derivative(unkdash(1,3),ff1,1)
      call derivative(unkdash(1,3),ff2,2)
      call derivative(unkdash(1,3),ff3,3)
      do ip=1,npoin
      ff1(ip)=ff1(ip)/vol(ip)
      ff2(ip)=ff2(ip)/vol(ip)
      ff3(ip)=ff3(ip)/vol(ip)
      enddo
      do ip=1,npoin
      uuu1(ip,3)=unkdash(ip,3)-timestp(ip)*
     &(uuu(ip,1)*ff1(ip)+uuu(ip,2)*ff2(ip)+uuu(ip,3)*ff3(ip))
      uuu1(ip,1)=uuu1(ip,1)/(gmm(ip)*cosa(ip))
      uuu1(ip,2)=uuu1(ip,2)/gmm(ip)
      pprpot(ip)=0.
        ff1(ip)=0.
        ff2(ip)=0.
        ff3(ip)=0.
      enddo
            call gcrk(pprpot,uuu1(1,1),uuu1(1,2),uuu1(1,3),
     &                           ff1(1),ff2(1),ff3(1),0)
      call prfr(pprpot,uuu1(1,1),uuu1(1,2),uuu1(1,3),
     &                   ff1(1),ff2(1),ff3(1),0)
      do ip=1,npoin
c  advective velocity at n+1/2
      uuu(ip,1)=(uuu(ip,1)+ff1(ip))*0.5*rh0(ip)
      uuu(ip,2)=(uuu(ip,2)+ff2(ip))*0.5*rh0(ip)
      uuu(ip,3)=(uuu(ip,3)+ff3(ip))*0.5*rh0(ip)
c  low order solution at n+1
      unkdash(ip,1)=ff1(ip)*gmm(ip)*cosa(ip)
      unkdash(ip,2)=ff2(ip)*gmm(ip)
      unkdash(ip,3)=ff3(ip)
      enddo
      return
      end
c----------------------

c
      subroutine thprof
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/ctherm/ rg,cp,cap,st,g,th00,tt00,pr00,rh00,u00,v00,u0z,v0z
      if(lipps.eq.1) then
      do 1 ip=1,npoin
    1 th0(ip)=th00*exp(st*coord(ip,3))
      else
      do 2 ip=1,npoin
    2 th0(ip)=th00
      endif
      return
      end
c
      subroutine rhprof
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      common/frflo/afree(5),p0inf,rh0inf,u0inf,v0inf,w0inf
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      common/profil/rh0(npoinx),th0(npoinx),the(npoinx),ue(npoinx),
     &ve(npoinx),lipps
      common/ctherm/ rgas,cp,cap,st,g,th00,
     &tt00,pr00,rh00,u00,v00,u0z,v0z
      rhmax=-1.e+06
      rhmin=1.e+06
      Rgas=287.04
      cp=3.5*Rgas
      cap=Rgas/cp
      capi=1./cap
      cs=g/(cp*tt00*st)
      rh0inf=1.
      sd=1.535e-04
      do 10 ip=1,npoin
      rh0(ip)=rh0inf
      if(lipps.eq.1)then
      exs=exp(-st*coord(ip,3))
      rh0(ip)=rh00*exs*(1.-cs*(1.-exs))**(capi-1.)
c      rh0(ip)=rh00*exp(-sd*coord(ip,3))
      endif
       rhmax=amax1(rhmax,rh0(ip))
      rhmin=amin1(rhmin,rh0(ip))
   10 continue
       write(6,*)'gestosc',rhmax,rhmin
      return
      end
c******************************************
c******************************************
      subroutine noise(frc)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
      dimension frc(npoinx,4)
      common/coin/ncoin,icoin(nfacesx*2,2),irecog1(nfacesx)
      common/n/npoin,nedge,nbface,nbpoin,nboun
      common/xyz/coord(npoinx,3)
      common/sf/isfera
      common/per/iper
      common/per2/  icoin2(nfacesx,2),ncoin2
      common/per4/  icoin4(nfacesx,4),ncoin4

c      randomf1(x,kk)=(x)
c      randomf2(x,kk)=(x-0.5)
c      randomf3(x,kk)=(x-0.5)*amax1(0.,1.-(kk-1)*dz/500.)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      inodes=1            ! different seeding on different processor
c     inodes=0            ! noise equivalent to single processor run
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Numerical Recipes in Fortran - Quick and Dirty Generators p.274-275
C       im         ia         ic           overflow
C      86436       1093      18254           2^27
C     117128       1277      24749           2^28
C     145800       3661      30809           2^29
C     139968       3877      29573           2^30
C     134456       8121      28411           2^31
C     233280       9301      49297           2^32
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      im=86436
      ia=1093
      ic=18254
        do ip=1,npoin

         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         frc(ip,1)=(randx-0.5)
     2                        *amax1(0.,1.-coord(ip,3)/500.)

         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         frc(ip,2)=(randx-0.5)


         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         frc(ip,3)=randx-0.5


        enddo
      if(isfera.eq.1)then
      do ic=1,ncoin
      i1=icoin(ic,1)
      i2=icoin(ic,2)
      p1=frc(i1,1)
      p2=frc(i2,1)
      frc(i1,1)=(p1+p2)*0.5
      frc(i2,1)=(p1+p2)*0.5
      p1=frc(i1,2)
      p2=frc(i2,2)
      frc(i1,2)=(p1+p2)*0.5
      frc(i2,2)=(p1+p2)*0.5
      p1=frc(i1,3)
      p2=frc(i2,3)
      frc(i1,3)=(p1+p2)*0.5
      frc(i2,3)=(p1+p2)*0.5
      enddo
      endif
      if(iper.eq.2)then
      do ic=1,ncoin
      i1=icoin(ic,1)
      i2=icoin(ic,2)
      p1=frc(i1,1)
      p2=frc(i2,1)
      frc(i1,1)=(p1+p2)*0.5
      frc(i2,1)=(p1+p2)*0.5
      p1=frc(i1,2)
      p2=frc(i2,2)
      frc(i1,2)=(p1+p2)*0.5
      frc(i2,2)=(p1+p2)*0.5
      p1=frc(i1,3)
      p2=frc(i2,3)
      frc(i1,3)=(p1+p2)*0.5
      frc(i2,3)=(p1+p2)*0.5
      enddo
      do ic=1,ncoin2
      i1=icoin2(ic,1)
      i2=icoin2(ic,2)
      p1=frc(i1,1)
      p2=frc(i2,1)
      frc(i1,1)=(p1+p2)*0.5
      frc(i2,1)=(p1+p2)*0.5
      p1=frc(i1,2)
      p2=frc(i2,2)
      frc(i1,2)=(p1+p2)*0.5
      frc(i2,2)=(p1+p2)*0.5
      p1=frc(i1,3)
      p2=frc(i2,3)
      frc(i1,3)=(p1+p2)*0.5
      frc(i2,3)=(p1+p2)*0.5
      enddo
      do ic=1,ncoin4
      i1=icoin4(ic,1)
      i2=icoin4(ic,2)
      i3=icoin4(ic,3)
      i4=icoin4(ic,4)
      p1=frc(i1,1)
      p2=frc(i2,1)
      p3=frc(i3,1)
      p4=frc(i4,1)
      frc(i1,1)=(p1+p2+p3+p4)*0.25
      frc(i2,1)=(p1+p2+p3+p4)*0.25
      frc(i3,1)=(p1+p2+p3+p4)*0.25
      frc(i4,1)=(p1+p2+p3+p4)*0.25
      p1=frc(i1,2)
      p2=frc(i2,2)
      p3=frc(i3,2)
      p4=frc(i4,2)
      frc(i1,2)=(p1+p2+p3+p4)*0.25
      frc(i2,2)=(p1+p2+p3+p4)*0.25
      frc(i3,2)=(p1+p2+p3+p4)*0.25
      frc(i4,2)=(p1+p2+p3+p4)*0.25
      p1=frc(i1,3)
      p2=frc(i2,3)
      p3=frc(i3,3)
      p4=frc(i4,3)
      frc(i1,3)=(p1+p2+p3+p4)*0.25
      frc(i2,3)=(p1+p2+p3+p4)*0.25
      frc(i3,3)=(p1+p2+p3+p4)*0.25
      frc(i4,3)=(p1+p2+p3+p4)*0.25
      enddo
      endif
      return
      end

        subroutine Restart_Out(u1,u2,u3,u4,u5,u6,p1,p2,p3,p4,
     &  n1,n2,ntime8)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (ncharlen=9)
        character(len=ncharlen) :: charI
        character(len=(ncharlen+12)) :: filename
        common/time/ntimet,nrstart,irestart,nre
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/Average_D/N_Av_flag,N_AvS,N_AvE,N_AvPre,N_av_flag_pre
        dimension u1(npoinx,3),u2(npoinx,3),u3(npoinx,3)
        dimension u4(npoinx,4),u5(npoinx,4),u6(npoinx,4),p1(npoinx)
        dimension p2(npoinx),p3(npoinx),p4(npoinx)



        nfid = 345
        write(charI,'(i9)')(ntime8+ntimet)
!        write(charI,lentgh)nfid
        filename ='Rstart_Out'//charI//'.d'
        open(nfid, file=filename)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do ip=1,npoin
          write(nfid,198)u1(ip,1),u1(ip,2),u1(ip,3)
          enddo
          do ip=1,npoin
          write(nfid,198)u2(ip,1),u2(ip,2),u2(ip,3)
          enddo
          do ip=1,npoin
           write(nfid,198)u3(ip,1), u3(ip,2),u3(ip,3)
          enddo
          do ip=1,npoin
          write(nfid,199)u4(ip,1),u4(ip,2),u4(ip,3),u4(ip,4)
          enddo
          do ip=1,npoin
          write(nfid,199)u5(ip,1),u5(ip,2),u5(ip,3),u5(ip,4)
          enddo
ccccccccccccccccccccccccccccc for average data cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc u6 and p4 cccccccccccccccccccccccccccc
            write(nfid,*)'average V'
          do ip=1,npoin
          write(nfid,199)u6(ip,1),u6(ip,2),u6(ip,3),u6(ip,4)
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do ip=1,npoin
          write(nfid,199)p1(ip),p2(ip),p3(ip),p4(ip)
          enddo
          write(nfid,*)n1,n2
          write(nfid,*)ntime8+ntimet
          write(nfid,*)N_AvS,N_AvE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
197   format(2f14.6)
198   format(3f14.6)
199   format(4f14.6)
         close(nfid)




        end



        subroutine Restart_Out_old(u1,u2,u3,u4,u5,p1,p2,p3,p4
     &       ,n1,n2,ntime8)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (ncharlen=9)
        character(len=ncharlen) :: charI
        character(len=(ncharlen+12)) :: filename
        common/time/ntimet,nrstart,irestart,nre
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/Average_D/N_Av_flag,N_AvS,N_AvE,N_AvPre,N_av_flag_pre
        dimension u1(npoinx,3),u2(npoinx,3),u3(npoinx,3)
        dimension u4(npoinx,4),u5(npoinx,4),p1(npoinx)
        dimension p2(npoinx),p3(npoinx),p4(npoinx)


        nfid = 345
        write(charI,'(i9)')(ntime8+ntimet)
!        write(charI,lentgh)nfid
        filename ='Rstart_Out'//charI//'.d'
        open(nfid, file=filename)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do ip=1,npoin
          write(nfid,198)u1(ip,1),u1(ip,2),u1(ip,3)
          enddo
          do ip=1,npoin
          write(nfid,198)u2(ip,1),u2(ip,2),u2(ip,3)
          enddo
          do ip=1,npoin
           write(nfid,198)u3(ip,1), u3(ip,2),u3(ip,3)
          enddo
          do ip=1,npoin
          write(nfid,199)u4(ip,1),u4(ip,2),u4(ip,3),u4(ip,4)
          enddo
          do ip=1,npoin
          write(nfid,199)u5(ip,1),u5(ip,2),u5(ip,3),u5(ip,4)
          enddo
          do ip=1,npoin
          write(nfid,199)p1(ip),p2(ip),p3(ip),p4(ip)
          enddo
          write(nfid,*)n1,n2
          write(nfid,*)ntime8+ntimet
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
197   format(2f14.6)
198   format(3f14.6)
199   format(4f14.6)
         close(nfid)




        end





      subroutine Restart_In_old(u1,u2,u3,u4,u5,p1,p2,p3,p4,n1,n2)

!    for the restart file that does not contain averaged data


      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (ncharlen=9)
        character(len=ncharlen) :: charI
!        character(len=1) :: Temp
!        character(len=4) :: lentgh
        character(len=(ncharlen+12)) :: filename
        common/time/ntimet,nrstart,irestart,nre
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/Average_D/N_Av_flag,N_AvS,N_AvE,N_AvPre,N_av_flag_pre
        dimension u1(npoinx,3),u2(npoinx,3),u3(npoinx,3)
        dimension u4(npoinx,4),u5(npoinx,4),p1(npoinx)
        dimension p2(npoinx),p3(npoinx),p4(npoinx)




        nfid = 346
        Rewind(346)
        write(charI,'(i9)')ntimet
        filename ='Rstart_Out'//charI//'.d'
        open(nfid, file=filename, access='sequential',status='old')


          do ip=1,npoin
          read(nfid,198)u1(ip,1),u1(ip,2),u1(ip,3)
          enddo
          do ip=1,npoin
          read(nfid,198)u2(ip,1),u2(ip,2),u2(ip,3)
          enddo
          do ip=1,npoin
           read(nfid,198)u3(ip,1), u3(ip,2),u3(ip,3)
          enddo
          do ip=1,npoin
          read(nfid,199)u4(ip,1),u4(ip,2),u4(ip,3),u4(ip,4)
          enddo
          do ip=1,npoin
          read(nfid,199)u5(ip,1),u5(ip,2),u5(ip,3),u5(ip,4)
          enddo
          do ip=1,npoin
          read(nfid,199)p1(ip),p2(ip),p3(ip),p4(ip)
          enddo
          read(nfid,*)n1,n2
c            endif
197   format(2f14.6)
198   format(3f14.6)
199   format(4f14.6)


        close(nfid)
        return

        end




      subroutine Restart_In(u1,u2,u3,u4,u5,u6,p1,p2,p3,p4,n1,n2)
      parameter (nedgex=4000000,npoinx=800000,nfacesx=110000)
        parameter (ncharlen=9)
        character(len=ncharlen) :: charI
!        character(len=1) :: Temp
!        character(len=4) :: lentgh
        character(len=(ncharlen+12)) :: filename
        common/time/ntimet,nrstart,irestart,nre
        common/n/npoin,nedge,nbface,nbpoin,nboun
        common/Average_D/N_Av_flag,N_AvS,N_AvE,N_AvPre,N_av_flag_pre
        dimension u1(npoinx,4),u2(npoinx,4),u3(npoinx,4)
        dimension u4(npoinx,4),u5(npoinx,4),u6(npoinx,4)
        dimension p1(npoinx),p2(npoinx),p3(npoinx),p4(npoinx)




        nfid = 346
        Rewind(346)
        write(charI,'(i9)')ntimet
        filename ='Rstart_Out'//charI//'.d'
        open(nfid, file=filename, access='sequential',status='old')


          do ip=1,npoin
          read(nfid,198)u1(ip,1),u1(ip,2),u1(ip,3),u1(ip,4)
          enddo
          do ip=1,npoin
          read(nfid,198)u2(ip,1),u2(ip,2),u2(ip,3),u2(ip,4)
          enddo
          do ip=1,npoin
           read(nfid,198)u3(ip,1), u3(ip,2),u3(ip,3),u3(ip,4)
          enddo
          do ip=1,npoin
          read(nfid,199)u4(ip,1),u4(ip,2),u4(ip,3),u4(ip,4)
          enddo
          do ip=1,npoin
          read(nfid,199)u5(ip,1),u5(ip,2),u5(ip,3),u5(ip,4)
          enddo

ccccccccccccccccccccccccccccc for average data cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc u6 and p4 cccccccccccccccccccccccccccc
          read(nfid,*)
          do ip=1,npoin
          read(nfid,199)u6(ip,1),u6(ip,2),u6(ip,3),u6(ip,4)
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do ip=1,npoin
          read(nfid,199)p1(ip),p2(ip),p3(ip),p4(ip)
          enddo
          read(nfid,*)n1,n2
          read(nfid,*)
!------   wish to restart and average from previous averaged result, N_av_flag_pre=1
!------   if there is no previous averaged result N_AvPre=0 else N_AvPre= the previous averaged iterations ---
!------   do not wish to restart from previous averaged result, N_av_flag_pre=0, can be set in inputre.d -----
          read(nfid,*)n_avstemp,n_avetemp
          if (N_av_flag_pre .eq. 1)then
            if(n_avetemp .ne. ntimet)then
            write(*,*)'the end of previous averaged iteration does',
     +    'not match the start of current run'
            stop
            endif
          if((n_avetemp-n_avstemp) .gt. 0)then
          N_AvPre=n_avetemp-n_avstemp
          else
          N_AvPre=0
          endif

          else
          do ip=1,npoin
            do i=1,4
              u6(ip,i)=0.
            enddo
            p4(ip)=0.
          enddo
          N_AvPre=0
          endif
c            endif
197   format(2f14.6)
198   format(3f14.6)
199   format(4f14.6)


        close(nfid)
        return

        end
