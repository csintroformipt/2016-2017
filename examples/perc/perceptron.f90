    program perceptron
    implicit none

    interface
        subroutine get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, vec, sig, avg, a)
            integer L, i
            real smin, smax,  avgmin, avgmax, amin, amax, vec(:), sig, avg, a
        end subroutine

        subroutine set_vec_gaussoid (L, vec, avg, sig, A)
            integer L
            real  vec(:), sig, avg, A
        end subroutine


        subroutine get_fit(vec, va, vr, wmean, wsig, wa , L, L1, L1t, Ncells, estmean, estsig, esta)
            integer L, Ncells
            real vec(:), va(:,:), vr, wmean(:), wsig(:), wa(:), L1(:), L1t(:)
            real estmean, estsig, esta
        end subroutine get_fit

        subroutine write_model (L, Ncells, va, vr, wmean, wsig,wa)
            integer L, Ncells
            real va(:,:), vr, wmean(:), wsig(:), wa(:)
        end subroutine write_model

        subroutine read_model (filename, L, Ncells, va, vr, wmean, wsig, wa)
            integer L, Ncells, va(:,:), vr, wmean(:), wsig(:), wa(:)
            integer un
            character filename(:)
        end subroutine read_model
    end interface


    integer, parameter:: L=100, Ncells=10, NLearningCicles=1000000
    integer,parameter:: Na=20, Nmean=20,Nsig=20
    integer,parameter:: NStudy=Na*Nmean*Nsig
    integer i, j, k, t
    real rnd
    real, allocatable:: study(:,:)
    real, allocatable:: stmean(:), stsig(:), sta(:)
    real, allocatable:: vec(:)
    real sig, avg, A, r
    real smin,smax, avgmin, avgmax, amax, amin      !tutor-set parameters
    real vr
    real, allocatable:: va(:,:)!, vr(:)       !a-layer weights and r-layer weights
    real, allocatable:: L1(:), L1t(:)        !first layer a-result and r-result
    real, allocatable:: wmean(:),wsig(:), wa(:)
    real estmean, estsig, esta
    real etha
    real ermean, ersig, era, ernorm
    real,allocatable:: erL1(:)

    real wampl1,wamplmean,wamplsig,wampla
    real, allocatable:: tvec(:)
    
    allocate (tvec(L))
    allocate(study(L,NStudy))
    allocate(stmean(NStudy))
    allocate (stsig(NStudy))
    allocate (sta(Nstudy))
    allocate (vec(L))

    allocate (va(L,Ncells))
    !allocate (vr(Ncells))

    allocate (L1(Ncells))
    allocate (L1t(Ncells))

    allocate (wmean(Ncells))
    allocate (wsig(Ncells))
    allocate (wa(Ncells))

    allocate (erL1(Ncells))

    smin = 20
    smax = 30
    avgmin = 20
    avgmax = 80
    amin = 1
    amax = 10

    ! fixed tutor set
    t=1
    do i=1, Na
        do j=1,Nmean
            do k=1,Nsig
                sig=smin+(smax-smin)*k/Nsig
                avg=avgmin+(avgmax-avgmin)*j/Nmean
                A=amin+(amax-amin)*i/Na
                call set_vec_gaussoid (L, tvec, avg, sig, A)
                study(:,t)=tvec
                stsig(t)=sig
                stmean(t)=avg
                sta(t)=a
                t=t+1
            end do
        end do
        write(*,*) 't=',t
    end do

    !!get randomized tutor set
    !!do j = 1, NStudy
    !do j=1,Na*Nmean*Nsig
    !    call get_gaussoid (smin, smax, L, tvec, sig, avg)
    !
    !    study(:,j)=tvec
    !    stsig(j)=sig
    !    stmean(j)=avg
    !    write(*,*) "avg=",avg
    !    write(*,*) "sig=",sig
    !    write(*,*)
    !end do



    !init weights
    !va=1.0
    wampl1=1.0
    wamplmean=1.0
    wamplsig=1.0
    wampla=1.0

    !set initial weights distribution amplitude [-wampl,wampl)
    do i=1,Ncells
        do j=1,L
            call RANDOM_NUMBER(rnd)
            va(j,i)=wampl1*2*(rnd-0.5)
        end do
        call RANDOM_NUMBER(rnd)
        wmean(i)=wamplmean*2*(rnd-0.5)
        call RANDOM_NUMBER(rnd)
        wsig(i)=wamplsig*2*(rnd-0.5)
        call RANDOM_NUMBER(rnd)
        wa(i)=wampla*2*(rnd-0.5)
    end do
    vr = 0.5
    etha = 1e-5

    ernorm = 0.0
    do t=1, NLearningCicles

        !backprop with mentor's set
        do k=1, Nstudy

            call get_fit(study(:,k), va, vr, wmean, wsig, wa, L, L1, L1t, Ncells, estmean, estsig, esta)

            !output layer errors
            ermean= (estmean-stmean(k))
            ersig = (estsig -stsig(k))
            era   = (esta   -sta(k))
            !ernorm = sqrt(ermean**2+ersig**2 )

            !modify weights
            !output layer errors are ermean and ersig
            !backpropagate errors
            do i=1,Ncells
                erL1(i)=2*vr*L1t(i)*(1-L1t(i))*(wmean(i)*ermean+wsig(i)*ersig+wa(i)*era)
            end do
            ! correct weights
            do i=1,Ncells
                !L1 weights
                do j =1,L
                    va(j,i)=va(j,i)-etha*study(j,k)*erL1(i)
                end do
                !output layer weights
                wmean(i)=wmean(i)-etha*L1t(i)*ermean
                wsig(i)=wsig(i)  -etha*L1t(i)*ersig
                wa(i)  =wa(i)    -etha*L1t(i)*era
            end do
        end do


        !control
        call get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, tvec, sig, avg, a)
        call get_fit(tvec, va, vr, wmean, wsig, wa, L, L1, L1t, Ncells, estmean, estsig, esta)
        
        ernorm = ernorm + ((estmean-avg)**2+(estsig-sig)**2+(esta-a)**2)
        if (mod(t,100)==0) then
        write(*,*) t
        write(*,*) "test_mean=",avg, "est_mean=" ,estmean, "mean_diff=", avg-estmean
        write(*,*) "test_sig =",sig ,"est_sig =" ,estsig,  "sig_diff =", sig-estsig
        write(*,*) "test_a   =",a   ,"est_a   =" ,esta,    "a_diff   =", a-esta
        write(*,*) "diff_norm_last100=",ernorm
        write(*,*)
        ernorm = 0
        end if
        if (mod(t,1000)==0) then
            etha=etha*0.99
        end if
    end do

    call write_model (L, Ncells, va, vr, wmean, wsig, wa)
    !!!!! end program !!!!!

    end program perceptron

    subroutine get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, vec, sig, avg, A)
        integer L
        real smin, smax,  avgmin, avgmax, amin, amax, vec(:), sig, avg, A
        
        integer i
        real r
        call RANDOM_NUMBER(r)
        sig = smin+(smax-smin)*r
        call RANDOM_NUMBER(r)
        avg = avgmin+(avgmax-avgmin)*r
        call RANDOM_NUMBER(r)
        A = amin + (amax-amin)*r
        do i = 1, L
            vec(i)=A*exp(-(avg-i)*(avg-i)*0.5/(sig*sig))
        end do
    end subroutine

    subroutine set_vec_gaussoid (L, vec, avg, sig, A)
        integer L
        real  vec(:), sig, avg, A
        
        integer i
        real r
        do i = 1, L
            vec(i)=A*exp(-(avg-i)*(avg-i)*0.5/(sig*sig))
        end do
    end subroutine


    subroutine get_fit(vec, va, vr, wmean, wsig, wa , L, L1, L1t, Ncells, estmean, estsig, esta)
    integer L, Ncells
    real vec(:), va(:,:), vr, wmean(:), wsig(:), wa(:), L1(:), L1t(:)
    real estmean, estsig, esta

    integer i, j

    ! get output
    !L0->L1
    do i=1,Ncells
        L1=0
        do j =1,L
            L1(i)=L1(i)+vec(j)*va(j,i)
        end do
        L1t(i)=1./(1+exp(-2*vr*L1(i)))
    end do

    !L1-> output
    estmean = 0.
    estsig  = 0.
    esta    = 0.
    do i=1,Ncells
        estmean=estmean+L1t(i)*wmean(i)
        estsig =estsig +L1t(i)*wsig(i)
        esta   =esta   +L1t(i)*wa(i)
    end do
    end subroutine get_fit

    subroutine write_model (L, Ncells, va, vr, wmean, wsig,wa)
    integer L, Ncells
    real va(:,:), vr, wmean(:), wsig(:), wa(:)
    integer un
    un=1
    open(un,file="latest_mode.dat",mode='write')
    write(un,*) L, Ncells
    write(un,*) va
    write(un,*) vr
    write(un,*) wmean
    write(un,*) wsig
    write(un,*) wa
    close(un)
    end subroutine write_model

    subroutine read_model (filename, L, Ncells, va, vr, wmean, wsig, wa)
    integer L, Ncells
    real va(:,:), vr, wmean(:), wsig(:), wa(:)
    integer un
    character filename(:)

    un=1
    open(un,file=filename,mode='read')
    write(un,*) L, Ncells
    write(un,*) va
    write(un,*) vr
    write(un,*) wmean
    write(un,*) wsig
    write(un,*) wa
    close(un)
    end subroutine read_model
