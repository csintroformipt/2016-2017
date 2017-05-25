    program perceptron
    implicit none

    interface
    subroutine get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, vec, sig, avg, a)
    integer L, i
    real smin, smax,  avgmin, avgmax, amin, amax, vec(:), sig, avg, a
    end subroutine

    subroutine get_fit(vec, va, vr, wmean, wsig, wa , L, L1, L1t, Ncells, estmean, estsig, esta)
    integer L, Ncells
    real vec(:), va(:,:), vr, wmean(:), wsig(:), wa(:), L1(:), L1t(:)
    real estmean, estsig, esta
    end subroutine get_fit

    subroutine read_model (filename, L, Ncells, va, vr, wmean, wsig, wa)
    integer L, Ncells
    real va(:,:), vr, wmean(:), wsig(:), wa(:)
    integer un
    character(255) filename
    end subroutine read_model
    
    subroutine read_model_size (filename, L, Ncells)
    integer L, Ncells
    integer un
    character(255) filename

    end subroutine read_model_size
    
    end interface


    integer L, Ncells

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

    character(255) fn
    fn="latest_model.dat"
    
    call read_model_size(fn,L,Ncells)
    
    allocate (tvec(L))

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

    vr = 0.5

    ernorm = 0.0
    
    call read_model (fn, L, Ncells, va, vr, wmean, wsig, wa)
    write(*,*) "readin is done"
    call get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, tvec, sig, avg, a)
    call get_gaussoid (smin, smax, avgmin, avgmax, amin, amax, L, tvec, sig, avg, a)
    call get_fit(tvec, va, vr, wmean, wsig, wa, L, L1, L1t, Ncells, estmean, estsig, esta)

    ernorm = ernorm + ((estmean-avg)**2+(estsig-sig)**2+(esta-a)**2)
        write(*,*) t
        write(*,*) "test_mean=",avg, "est_mean=" ,estmean, "mean_diff=", avg-estmean
        write(*,*) "test_sig =",sig ,"est_sig =" ,estsig,  "sig_diff =", sig-estsig
        write(*,*) "test_a   =",a   ,"est_a   =" ,esta,    "a_diff   =", a-esta
        write(*,*) "diff_norm_last100=",ernorm
        write(*,*)
        ernorm = 0



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

    
    subroutine read_model_size (filename, L, Ncells)
    integer L, Ncells
    integer un
    character(255) filename

    un=1
    open(un,file="latest_model.dat")
    read(un,*) L, Ncells
    close(un)
    end subroutine read_model_size
    
    subroutine read_model (filename, L, Ncells, va, vr, wmean, wsig, wa)
    integer L, Ncells
    real va(:,:), vr, wmean(:), wsig(:), wa(:)
    integer un
    character(255) filename

    un=1
    open(un,file="latest_model.dat")
    write(*,*) 123
    read(un,*) L, Ncells
    write(*,*) 123
    read(un,*) va
    read(un,*) vr
    read(un,*) wmean
    read(un,*) wsig
    read(un,*) wa
    close(un)
    end subroutine read_model
