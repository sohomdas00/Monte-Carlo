! program to simulate Ising Model using Metropolis Monte Carlo
     module parameters
       implicit none
       integer*4,parameter :: l=64        ! system size
       integer*4,dimension(l,l) :: s
       integer*4 :: seed1,seed2,seed3
       real*8 :: temp,en,mag
       integer*4 :: icycl
     end module parameters

    program mc
    use parameters
    implicit none
    integer*4,parameter::ncycl=10000
    
    open(25,file='seed.dat',status='old')
    read(25,*) seed1
    read(25,*) seed2
    read(25,*) seed3

    open(26,file='temp.dat',status='old')
    read(26,*) temp

    open(27,file='energy.dat',status='unknown')
    open(28,file='mag.dat',status='unknown')

    call init           ! initialises the lattice
    call energy
    call magnetisation
    write(27,*) 0,en
    write(28,*) 0,mag

    DO icycl=1,ncycl
     call mcmove
     call energy
     call magnetisation 
     !if (mod(icycl,10).eq.0) call save_snap 
     write(27,*)icycl, en
     write(28,*)icycl, mag
    ENDDO
    END
!------------------------------------------------------------
    SUBROUTINE init
    use parameters
    implicit none
    integer*4 :: i,j
    real*8 :: ran3,r
    
    open(56,file='snaps/snap_0000',status='unknown')

    do i=1,l
    do j=1,l
      r = ran3(seed3)
      if (r.lt.0.5) then
         s(i,j)=-1
      else
         s(i,j)=+1
      endif
      write(56,*) i,j,s(i,j)
    enddo
    enddo

   return
   end
       
!------------------------------------------------------------
    SUBROUTINE mcmove
    use parameters
    implicit none
    integer*4 :: i,j,ox,oy
    integer*4 :: oxr,oxl,oyr,oyl
    real*8 :: ran3,nc,engy,r
    
    do i=1,l
      do j=1,l
        ox = int(ran3(seed1)*real(l)) + 1
        oy = int(ran3(seed2)*real(l)) + 1
        if (ox.gt.l)ox=1
        if (oy.gt.l)oy=1
        oxr=ox+1
        oxl=ox-1
        oyr=oy+1
        oyl=oy-1
        if (oxr.gt.l)oxr=1
        if (oyr.gt.l)oyr=1
        if (oxl.lt.1)oxl=l
        if (oyl.lt.1)oyl=1
        nc = real(s(oxr,oy)+s(ox,oyr)+s(oxl,oy)+s(ox,oyl))
        engy = 2.0*real(s(ox,oy))*nc
        r = ran3(seed3)
        if (r.lt.exp(-engy/temp)) then
          s(ox,oy) = -s(ox,oy)  ! flip spin
        endif
      enddo
    enddo
    return 
    end
!------------------------------------------------------------
    SUBROUTINE save_snap
    use parameters
    implicit none
    integer*4 :: i,j
    character(len=100) :: name1
    
    write(name1,'(a,i0.4)') 'snaps/snap_',icycl
    open(34,file=name1,status='unknown')

    do i=1,l
    do j=1,l
      write(34,*) i, j, s(i,j)
    enddo
    enddo

    return
    end
!----------------------------------------------------------
    SUBROUTINE energy
    use parameters
    implicit none
    integer*4 :: i,j
    integer*4 :: ir,il,jr,jl
    real*8 :: enn

    en = 0.0
    do i=1,l
    do j=1,l
      ir=i+1
      il=i-1
      jr=j+1
      jl=j-1
      if (ir.gt.l)ir=1
      if (jr.gt.l)jr=1
      if (il.lt.1)il=l
      if (jl.lt.1)jl=l
      enn = real(s(ir,j)+s(il,j)+s(i,jr)+s(i,jl))
      enn = -enn*real(s(i,j))
      en = en + enn
    enddo
    enddo
    en = en/(real(l*l)*4.0) 
    return 
    end
!-----------------------------------------------------------
    SUBROUTINE magnetisation
    use parameters
    implicit none
    integer*4 :: i,j
    
    mag=0.0

    do i=1,l
    do j=1,l
      mag=mag+real(s(i,j))
    enddo
    enddo
    
    mag = mag/real(l*l)

    return 
    end
!-----------------------------------------------------------      
        FUNCTION ran3(idum)
        INTEGER*4 idum
        INTEGER*4 MBIG,MSEED,MZ
  !      REAL MBIG,MSEED,MZ
        REAL*8 ran3,FAC
        PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  !     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
        INTEGER*4 i,iff,ii,inext,inextp,k
        INTEGER*4 mj,mk,ma(55)
  !     REAL mj,mk,ma(55)
        SAVE iff,inext,inextp,ma
        DATA iff /0/
        if(idum.lt.0.or.iff.eq.0)then
          iff=1
          mj=MSEED-iabs(idum)
          mj=mod(mj,MBIG)
          ma(55)=mj
          mk=1
          do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
  11      continue
          do 13 k=1,4
            do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
  12        continue
  13      continue
          inext=0
          inextp=31
          idum=1
        endif
        inext=inext+1
        if(inext.eq.56)inext=1
        inextp=inextp+1
        if(inextp.eq.56)inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj.lt.MZ)mj=mj+MBIG
        ma(inext)=mj
        ran3=mj*FAC
        return
        END

