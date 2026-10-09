c----------------------------------------------------------------
c
c Converts protein primary sequence to 3-letter code where the
c labels F, D, and P predict protein regions that are folded (F),
c intrinsically disordered (D), or phase-separating intrinsically
c disordered (P), respectively. Also the summed P classified distance
c is calculated with and without Uπ (π–π and π-cation contributions)
c and Uq (charge-based contributions).
c
c Program input:
c
c primary sequence, no gaps, restricted to the 20 common amino acid types
c minimum length = 25
c maximum length = 10000
c
c Program output:
c
c Converted sequence
c Summed classifier distance of P-labeled windows
c Summed classifier distance of P-labeled windows with Uπ and Uq corrections
c CSV file of Residue number, Amino Acid type, Residue label, Classifier Distance, Residue label (w/ Uπ Uq extension), Classifier Distance (w/ Uπ Uq extension)
c
c To use Parse_v2:
c
c Compile Parse_v2.f using a fortran compiler. The program was tested using GNU Fortran (Homebrew GCC 11.3.0_2) 11.3.0.
c Run the program at a command prompt, by using the executable followed by the protein primary sequence.
c For example "./a.out SEQUENCESEQUENCESEQUENCESEQUENCE"
c
c
c STW 03/30/2022
c
c----------------------------------------------------------------

      program Parse_v2
      implicit none

      integer length,num,npep,window_size,middle_position,
     &        num_ala_w,num_cys_w,num_asp_w,num_glu_w,num_phe_w,
     &        num_gly_w,num_his_w,num_ile_w,num_lys_w,num_leu_w,
     &        num_met_w,num_asn_w,num_pro_w,num_gln_w,num_arg_w,
     &        num_ser_w,num_thr_w,num_val_w,num_trp_w,num_tyr_w,
     &        count_regions,count_p,count_w,region,
     &        p_region_start,p_region_end,p_start(10000),
     &        p_end(10000),count_regions_p,count_regions_d,
     &        count_regions_f,d_start(10000),d_end(10000),
     &        f_start(10000),f_end(10000),count_domains,low,
     &        domain(10000),i,j,jj,jjj,jjjj
      character code_inp*11000,code(11000)*3,classification(11000)*3,
     &          classification_pi_q(11000)*3,
     &          classification_pi_q_csat(11000)*3
      real pppiia_h,pppiip_h,pppiig_h,pppiic_h,pppiid_h,pppiie_h,
     &     pppiif_h,pppiih_h,pppiii_h,pppiik_h,pppiil_h,pppiim_h,
     &     pppiin_h,pppiiq_h,pppiir_h,pppiis_h,pppiit_h,pppiiv_h,
     &     pppiiw_h,pppiiy_h,sum_ppii,v_exponent,fppii,nu_model,
     &     helix_a,helix_c,helix_d,helix_e,helix_f,helix_g,helix_h,
     &     helix_i,helix_k,helix_l,helix_m,helix_n,helix_p,helix_q,
     &     helix_r,helix_s,helix_t,helix_v,helix_w,helix_y,
     &     hydr_a,hydr_c,hydr_d,hydr_e,hydr_f,hydr_g,hydr_h,
     &     hydr_i,hydr_k,hydr_l,hydr_m,hydr_n,hydr_p,hydr_q,
     &     hydr_r,hydr_s,hydr_t,hydr_v,hydr_w,hydr_y,m,b,x,y,
     &     percent_p,percent_cutoff,net_charge_w,dist_norm(10000),
     &     rh_w,helix,hydr,ID_dist,PS_dist,p_dist_sum,U_pi,U_q,
     &     RY,RF,KY,KF,FY,lambda_1,lambda_2,ncpr_w,scd_w,
     &     p_pi_q_dist_sum,dist_norm_pi_q(10000),U_pi_csat,
     &     U_q_csat,p_pi_q_csat_dist_sum,dist_norm_pi_q_csat(10000)

c PPII bias measured in peptides by Hilser group; used to calculate Rh and then nu_model.
c Prot Sci 2013, vol 22, pgs 405-417, in a table in supplementary information

      pppiia_h=0.37
      pppiic_h=0.25
      pppiid_h=0.30
      pppiie_h=0.42
      pppiif_h=0.17
      pppiig_h=0.13
      pppiih_h=0.20
      pppiii_h=0.39
      pppiik_h=0.56
      pppiil_h=0.24
      pppiim_h=0.36
      pppiin_h=0.27
      pppiip_h=1.00
      pppiiq_h=0.53
      pppiir_h=0.38
      pppiis_h=0.24
      pppiit_h=0.32
      pppiiv_h=0.39
      pppiiw_h=0.25
      pppiiy_h=0.25

c Normalized frequency of alpha-helix (Tanaka-Scheraga, Macromolecules 10, 9-20 (1977); recalculated by Kidera as normalized frequencies).

      helix_a=1.42
      helix_c=0.73
      helix_d=1.01
      helix_e=1.63
      helix_f=1.16
      helix_g=0.50
      helix_h=1.20
      helix_i=1.12
      helix_k=1.24
      helix_l=1.29
      helix_m=1.21
      helix_n=0.71
      helix_p=0.65
      helix_q=1.02
      helix_r=1.06
      helix_s=0.71
      helix_t=0.78
      helix_v=0.99
      helix_w=1.05
      helix_y=0.67

c Structure-based interactivity scale used to calculate the hydrophobicity profile of a protein
c from its primary sequence. This intertactivity scale was determined from the residue contact matrix
c of single-domain globular proteins (Bastolla et al., Proteins 58, 22-30 (2005)).

      hydr_a=0.0728
      hydr_c=0.3557
      hydr_d=-0.0552
      hydr_e=-0.0295
      hydr_f=0.4201
      hydr_g=-0.0589
      hydr_h=0.0874
      hydr_i=0.3805
      hydr_k=-0.0053
      hydr_l=0.3819
      hydr_m=0.1613
      hydr_n=-0.0390
      hydr_p=-0.0492
      hydr_q=0.0126
      hydr_r=0.0394
      hydr_s=-0.0282
      hydr_t=0.0239
      hydr_v=0.2947
      hydr_w=0.4114
      hydr_y=0.3113


c Read input protein sequence.

      call get_command_argument(1,code_inp)
      if (len_trim(code_inp) == 0) then
      write(*,*)'no input argument, exiting program'
      stop
      endif

      length=len(code_inp)

c  convert any lower case letter to upper case.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

c  Determine sequence length.

      j=0
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 1
         j=j+1
         code(j)=code_inp(i:i)
1     continue
      enddo 
      npep=j

c restrict protein sequence to the 20 common amino acid types; check for numbers too

      do i=1,npep
      if (code(i).eq.'B'.or.code(i).eq.'J'.or.code(i).eq.'O'.or.
     &    code(i).eq.'U'.or.code(i).eq.'X'.or.code(i).eq.'Z'.or.
     &    code(i).eq.'1'.or.code(i).eq.'2'.or.code(i).eq.'3'.or.
     &    code(i).eq.'4'.or.code(i).eq.'5'.or.code(i).eq.'6'.or.
     &    code(i).eq.'7'.or.code(i).eq.'8'.or.code(i).eq.'9'.or.
     &    code(i).eq.'0') then
      write(*,*)' '
      write(*,*)'could not parse sequence'
      write(*,*)'sequence contains noncommon amino acid type'
      write(*,*)' '
      stop
      endif
      enddo

c Define window size.

      window_size=25

c if protein sequence is less than the window size, stop

      if (npep.lt.window_size) then
      write(*,*)' '
      write(*,*)'could not parse sequence'
      write(*,*)'input sequence is too short'
      write(*,*)'minimum sequence length is ',window_size
      write(*,*)' '
      stop
      endif

c if protein sequence is too long, stop

      if (npep.gt.10000) then
      write(*,*)' '
      write(*,*)'could not parse sequence'
      write(*,*)'input sequence is too long'
      write(*,*)'maximum sequence length is 10000'
      write(*,*)' '
      stop
      endif

      p_dist_sum=0.0
      p_pi_q_dist_sum=0.0
      p_pi_q_csat_dist_sum=0.0


c calculate hydrophobicity, nu_model, helix propensity, U_pi and Uq for each 25-residue window

      DO J=1,NPEP
         if (j.le.(NPEP-window_size+1)) then
         middle_position=j+window_size/2

      num_ala_w=0
      num_cys_w=0
      num_asp_w=0
      num_glu_w=0
      num_phe_w=0
      num_gly_w=0
      num_his_w=0
      num_ile_w=0
      num_lys_w=0
      num_leu_w=0
      num_met_w=0
      num_asn_w=0
      num_pro_w=0
      num_gln_w=0
      num_arg_w=0
      num_ser_w=0
      num_thr_w=0
      num_val_w=0
      num_trp_w=0
      num_tyr_w=0

      DO JJJ=J,J+window_size-1
         IF (CODE(JJJ).EQ.'A') THEN
            num_ala_w=num_ala_w+1
         endif
         IF (CODE(JJJ).EQ.'C') THEN
            num_cys_w=num_cys_w+1
         endif
         IF (CODE(JJJ).EQ.'D') THEN
            num_asp_w=num_asp_w+1
         endif
         IF (CODE(JJJ).EQ.'E') THEN
            num_glu_w=num_glu_w+1
         endif
         IF (CODE(JJJ).EQ.'F') THEN
            num_phe_w=num_phe_w+1
         endif
         IF (CODE(JJJ).EQ.'G') THEN
            num_gly_w=num_gly_w+1
         endif
         IF (CODE(JJJ).EQ.'H') THEN
            num_his_w=num_his_w+1
         endif
         IF (CODE(JJJ).EQ.'I') THEN
            num_ile_w=num_ile_w+1
         endif
         IF (CODE(JJJ).EQ.'K') THEN
            num_lys_w=num_lys_w+1
         endif
         IF (CODE(JJJ).EQ.'L') THEN
            num_leu_w=num_leu_w+1
         endif
         IF (CODE(JJJ).EQ.'M') THEN
            num_met_w=num_met_w+1
         endif
         IF (CODE(JJJ).EQ.'N') THEN
            num_asn_w=num_asn_w+1
         endif
         IF (CODE(JJJ).EQ.'P') THEN
            num_pro_w=num_pro_w+1
         endif
         IF (CODE(JJJ).EQ.'Q') THEN
            num_gln_w=num_gln_w+1
         endif
         IF (CODE(JJJ).EQ.'R') THEN
            num_arg_w=num_arg_w+1
         endif
         IF (CODE(JJJ).EQ.'S') THEN
            num_ser_w=num_ser_w+1
         endif
         IF (CODE(JJJ).EQ.'T') THEN
            num_thr_w=num_thr_w+1
         endif
         IF (CODE(JJJ).EQ.'V') THEN
            num_val_w=num_val_w+1
         endif
         IF (CODE(JJJ).EQ.'W') THEN
            num_trp_w=num_trp_w+1
         endif
         IF (CODE(JJJ).EQ.'Y') THEN
            num_tyr_w=num_tyr_w+1
         endif
      enddo

c calculate hydrophobicity

      hydr=0.0

      hydr=num_ala_w*hydr_a+num_cys_w*hydr_c
     &        +num_asp_w*hydr_d+num_glu_w*hydr_e
     &        +num_phe_w*hydr_f+num_gly_w*hydr_g
     &        +num_his_w*hydr_h+num_ile_w*hydr_i
     &        +num_lys_w*hydr_k+num_leu_w*hydr_l
     &        +num_met_w*hydr_m+num_asn_w*hydr_n
     &        +num_pro_w*hydr_p+num_gln_w*hydr_q
     &        +num_arg_w*hydr_r+num_ser_w*hydr_s
     &        +num_thr_w*hydr_t+num_val_w*hydr_v
     &        +num_trp_w*hydr_w+num_tyr_w*hydr_y
      hydr=hydr/real(window_size)

c Sector boundaries were defined by the mean and standard deviation in nu_model, helix
c propensity, and hydrophobicity calculated for each of the PS ID, ID, and
c folded sets. For the folded set, mean hydr is 0.1163898 ± 0.01679414
c
c Windows with hydr value less than 0.08280152 (mean - 2*sd) are classified as disordered (D or P).
c Windows with hydr value greater than or equal to 0.08280152 are classified as folded (F).
c
c dist_norm is the distance from the sector boundary, normalized by the distance to the mean.

      if(hydr.ge.(0.08280152)) then
        classification(middle_position)='F'
        classification_pi_q(middle_position)='F'
        classification_pi_q_csat(middle_position)='F'
        dist_norm(middle_position)=(hydr-0.08280152)/(0.01679414*2.0)
        dist_norm_pi_q(middle_position)=dist_norm(middle_position)
        dist_norm_pi_q_csat(middle_position)=dist_norm(middle_position)
        goto 10
      endif
 
c Calculate U_pi (trained against ∆h° data) that models π-π and cation-π effects
c Calculate U_pi_csat (trained against csat data)

         RY=0.0
         RF=0.0
         KY=0.0
         KF=0.0
         FY=0.0

         if (num_arg_w.ne.num_tyr_w) then
         RY=real(num_arg_w*num_tyr_w)/abs(real(num_arg_w-num_tyr_w))
         else
         RY=real(num_arg_w*num_tyr_w)
         endif
         if (num_arg_w.ne.num_phe_w) then
         RF=real(num_arg_w*num_phe_w)/abs(real(num_arg_w-num_phe_w))
         else
         RF=real(num_arg_w*num_phe_w)
         endif
         if (num_lys_w.ne.num_tyr_w) then
         KY=real(num_lys_w*num_tyr_w)/abs(real(num_lys_w-num_tyr_w))
         else
         KY=real(num_lys_w*num_tyr_w)
         endif
         if (num_lys_w.ne.num_phe_w) then
         KF=real(num_lys_w*num_phe_w)/abs(real(num_lys_w-num_phe_w))
         else
         KF=real(num_lys_w*num_phe_w)
         endif
         if (num_phe_w.ne.num_tyr_w) then
         FY=real(num_phe_w*num_tyr_w)/abs(real(num_phe_w-num_tyr_w))
         else
         FY=real(num_phe_w*num_tyr_w)
         endif

         U_pi=3.0*RY + 2.0*KY + 2.0*RF + 1.0*KF + 1.0*FY
         U_pi_csat=0.28*U_pi
         U_pi=0.137*U_pi

c Calculate U_q (trained against ∆h° data) that models charge effects by NCPR and SCD
c Calculate U_q_csat (trained against csat data)

      scd_w=0.0
      do jjj=j,j+window_size-1
         do jjjj=j,j+window_size-1
            if(jjjj.gt.jjj) then
            if(code(jjj).eq.'A')lambda_1=0.0
            if(code(jjj).eq.'C')lambda_1=0.0
            if(code(jjj).eq.'D')lambda_1=-1.0
            if(code(jjj).eq.'E')lambda_1=-1.0
            if(code(jjj).eq.'F')lambda_1=0.0
            if(code(jjj).eq.'G')lambda_1=0.0
            if(code(jjj).eq.'H')lambda_1=0.0
            if(code(jjj).eq.'I')lambda_1=0.0
            if(code(jjj).eq.'K')lambda_1=1.0
            if(code(jjj).eq.'L')lambda_1=0.0
            if(code(jjj).eq.'M')lambda_1=0.0
            if(code(jjj).eq.'N')lambda_1=0.0
            if(code(jjj).eq.'P')lambda_1=0.0
            if(code(jjj).eq.'Q')lambda_1=0.0
            if(code(jjj).eq.'R')lambda_1=1.0
            if(code(jjj).eq.'S')lambda_1=0.0
            if(code(jjj).eq.'T')lambda_1=0.0
            if(code(jjj).eq.'V')lambda_1=0.0
            if(code(jjj).eq.'W')lambda_1=0.0
            if(code(jjj).eq.'Y')lambda_1=0.0
            if(code(jjjj).eq.'A')lambda_2=0.0
            if(code(jjjj).eq.'C')lambda_2=0.0
            if(code(jjjj).eq.'D')lambda_2=-1.0
            if(code(jjjj).eq.'E')lambda_2=-1.0
            if(code(jjjj).eq.'F')lambda_2=0.0
            if(code(jjjj).eq.'G')lambda_2=0.0
            if(code(jjjj).eq.'H')lambda_2=0.0
            if(code(jjjj).eq.'I')lambda_2=0.0
            if(code(jjjj).eq.'K')lambda_2=1.0
            if(code(jjjj).eq.'L')lambda_2=0.0
            if(code(jjjj).eq.'M')lambda_2=0.0
            if(code(jjjj).eq.'N')lambda_2=0.0
            if(code(jjjj).eq.'P')lambda_2=0.0
            if(code(jjjj).eq.'Q')lambda_2=0.0
            if(code(jjjj).eq.'R')lambda_2=1.0
            if(code(jjjj).eq.'S')lambda_2=0.0
            if(code(jjjj).eq.'T')lambda_2=0.0
            if(code(jjjj).eq.'V')lambda_2=0.0
            if(code(jjjj).eq.'W')lambda_2=0.0
            if(code(jjjj).eq.'Y')lambda_2=0.0
            scd_w=scd_w+(lambda_1*lambda_2)*sqrt(real(abs(jjjj-jjj)))
            endif
         enddo
      enddo
      scd_w=scd_w/real(window_size)

      net_charge_w=abs(num_asp_w+num_glu_w
     &                  -num_lys_w-num_arg_w)

      ncpr_w=real(net_charge_w)/real(window_size)

      U_q=8.4*scd_w+5.6*ncpr_w
      U_q_csat=-16.0*scd_w+33.0*ncpr_w

c calculate nu_model

      sum_ppii=0.0

      sum_ppii=num_ala_w*pppiia_h+num_cys_w*pppiic_h
     &        +num_asp_w*pppiid_h+num_glu_w*pppiie_h
     &        +num_phe_w*pppiif_h+num_gly_w*pppiig_h
     &        +num_his_w*pppiih_h+num_ile_w*pppiii_h
     &        +num_lys_w*pppiik_h+num_leu_w*pppiil_h
     &        +num_met_w*pppiim_h+num_asn_w*pppiin_h
     &        +num_pro_w*pppiip_h+num_gln_w*pppiiq_h
     &        +num_arg_w*pppiir_h+num_ser_w*pppiis_h
     &        +num_thr_w*pppiit_h+num_val_w*pppiiv_h
     &        +num_trp_w*pppiiw_h+num_tyr_w*pppiiy_h

      fppii=sum_ppii/real(window_size)

c if the sequence is polyproline, then fppii = 1.0, which causes an error from log (1-fppii)
c to avoid this rare error
      if(fppii.eq.(1.0)) fppii=0.98

      v_exponent=0.503-0.11*log(1.0-fppii)

c multiplier of 4 on sequence length, putting nu_model into length-independent range (see Paiz et al 2021 JBC 297(5) 101343)

      rh_w=2.16*(real(4*window_size)**(v_exponent))
     &    +0.26*real(4*net_charge_w)
     &    -0.29*(real(4*window_size)**(0.5))
      nu_model=log(rh_w/2.16)/log(real(4*window_size))

c calculate helix propensity

      helix=0.0

      helix=num_ala_w*helix_a+num_cys_w*helix_c
     &        +num_asp_w*helix_d+num_glu_w*helix_e
     &        +num_phe_w*helix_f+num_gly_w*helix_g
     &        +num_his_w*helix_h+num_ile_w*helix_i
     &        +num_lys_w*helix_k+num_leu_w*helix_l
     &        +num_met_w*helix_m+num_asn_w*helix_n
     &        +num_pro_w*helix_p+num_gln_w*helix_q
     &        +num_arg_w*helix_r+num_ser_w*helix_s
     &        +num_thr_w*helix_t+num_val_w*helix_v
     &        +num_trp_w*helix_w+num_tyr_w*helix_y
      helix=helix/real(window_size)


c For the PS ID set, mean helix propensity is 0.9327272 ± 0.090295
c For the PS ID set, mean nu_model is 0.5416 ± 0.01997657
c For the ID set, mean helix propensity is 1.022552 ± 0.08171361
c For the ID set, mean nu_model is 0.5582901 ± 0.02200711 
c
c Boundary between P and D sectors is defined by the line
c
c y=(-0.244078945)*x + 0.7885823
c
c determined from the PS ID and ID set means and standard deviations, and 
c designed to bisect the overlapping set distributions.
c
c In a plot of helix propensity (x) versus nu_model (y), a window localized
c to the right of this boundary line is in the D sector. A window localized
c to the left of the boundary is in the P sector.

c distance from the boundary line to the ID set means
      m=-1.0/(-0.244078945)
      b=0.5582901-m*1.022552
      x=(b-0.7885823)/(-0.244078945-m)
      y=m*x+b
      ID_dist=sqrt((1.022552-x)*(1.022552-x)+
     &              (0.5582901-y)*(0.5582901-y))

c distance from the boundary line to the PS ID set means
      m=-1.0/(-0.244078945)
      b=0.5416-m*0.9327272
      x=(b-0.7885823)/(-0.244078945-m)
      y=m*x+b
      PS_dist=sqrt((0.9327272-x)*(0.9327272-x)+
     &              (0.5416-y)*(0.5416-y))

c Below, x and y define the point on the boundary line that makes a
c perpendicular when paired with the point defined by the window values
c of nu_model and helix propensity.
c Below, m and b define the equation of this perpendicular line. 
c Because the line is perpendicular to the boundary, its slope will be
c the negative reciprocal of the boundary slope.
      m=-1.0/(-0.244078945)
      b=nu_model-m*helix
c x and y define the intersect point of the boundary (y=(-0.244078945)*x + 0.7885823)
c and perpendicular (y=mx+b). Two equations, two unknowns, thus easily solved.
      x=(b-0.7885823)/(-0.244078945-m)
      y=m*x+b

      if(((nu_model-0.7885823)/(-0.244078945)).le.(helix)) then
         classification(middle_position)='D'
         classification_pi_q(middle_position)='D'
         classification_pi_q_csat(middle_position)='D'
         dist_norm(middle_position)=
     &      sqrt((helix-x)*(helix-x)+(nu_model-y)*(nu_model-y))/
     &      ID_dist
         dist_norm_pi_q(middle_position)=dist_norm(middle_position)
         dist_norm_pi_q_csat(middle_position)=
     &      dist_norm(middle_position)

         if (dist_norm(middle_position).lt.(U_pi+U_q)) then
            classification_pi_q(middle_position)='P'
            p_pi_q_dist_sum=p_pi_q_dist_sum
     &      +U_pi+U_q-dist_norm(middle_position)
         dist_norm_pi_q(middle_position)=U_pi+U_q
     &                 -dist_norm(middle_position)
         endif

         if (dist_norm(middle_position).lt.(U_pi_csat+U_q_csat)) then
            classification_pi_q_csat(middle_position)='P'
            p_pi_q_csat_dist_sum=p_pi_q_csat_dist_sum
     &      +U_pi_csat+U_q_csat-dist_norm(middle_position)
         dist_norm_pi_q_csat(middle_position)=U_pi_csat+U_q_csat
     &                 -dist_norm(middle_position)
         endif

         goto 10
      else
         classification(middle_position)='P'
         classification_pi_q(middle_position)='P'
         classification_pi_q_csat(middle_position)='P'
         dist_norm(middle_position)=
     &      sqrt((helix-x)*(helix-x)+(nu_model-y)*(nu_model-y))/
     &      PS_dist
         dist_norm_pi_q(middle_position)=U_pi+U_q
     &                 +dist_norm(middle_position)
         dist_norm_pi_q_csat(middle_position)=U_pi_csat+U_q_csat
     &                 +dist_norm(middle_position)
            p_dist_sum=p_dist_sum+dist_norm(middle_position)
            p_pi_q_dist_sum=p_pi_q_dist_sum
     &                 +dist_norm_pi_q(middle_position)
            p_pi_q_csat_dist_sum=p_pi_q_csat_dist_sum
     &                 +dist_norm_pi_q_csat(middle_position)
         goto 10
      endif

10    continue
      endif
      enddo

      do j=1,window_size/2
         classification(j)=classification((window_size/2)+1)
         classification_pi_q(j)=classification_pi_q((window_size/2)+1)
         classification_pi_q_csat(j)=
     &      classification_pi_q_csat((window_size/2)+1)
         dist_norm(j)=dist_norm((window_size/2)+1)
         dist_norm_pi_q(j)=dist_norm_pi_q((window_size/2)+1)
         dist_norm_pi_q_csat(j)=dist_norm_pi_q_csat((window_size/2)+1)
      enddo
      do j=npep-(window_size/2)+1,npep
         classification(j)=classification(npep-(window_size/2))
         classification_pi_q(j)=
     &      classification_pi_q(npep-(window_size/2))
         classification_pi_q_csat(j)=
     &      classification_pi_q_csat(npep-(window_size/2))
         dist_norm(j)=dist_norm(npep-(window_size/2))
         dist_norm_pi_q(j)=dist_norm_pi_q(npep-(window_size/2))
         dist_norm_pi_q_csat(j)=
     &      dist_norm_pi_q_csat(npep-(window_size/2))
      enddo


      write(*,*)' '
      write(*,'(11000a1)')(classification(j),j=1,npep)

      write(*,*)' '
      write(*,'("Including Uπ + Uq extension (∆h° trained):")')
      write(*,'(11000a1)')(classification_pi_q(j),j=1,npep)

      write(*,*)' '
      write(*,'("Including Uπ + Uq extension (csat trained):")')
      write(*,'(11000a1)')(classification_pi_q_csat(j),j=1,npep)

      write(*,*)' '
      write(*,'("Sequence length ",i6)')npep
      write(*,'("∑ classifier distance of P-labeled windows ",f10.3)')
     & p_dist_sum
      write(*,'("∑ classifier distance of P-labeled windows +Uπ +Uq ",
     &      "(∆h° trained)",f10.3)')p_pi_q_dist_sum
      write(*,'("∑ classifier distance of P-labeled windows +Uπ +Uq ",
     &      "(csat trained)",f10.3)')p_pi_q_csat_dist_sum

      open (7,file='residue_label_classifier_distance.csv')
      do j=1,npep
      write(7,'(i6,", ",a1,", ",a1,",",f8.3,", ",a1,",",f8.3,
     & ", ",a1,",",f8.3)')j,
     & code(j),classification(j),dist_norm(j),
     & classification_pi_q(j),dist_norm_pi_q(j),
     & classification_pi_q_csat(j),dist_norm_pi_q_csat(j)
      enddo
      close(7)

      percent_cutoff=0.90
c find PS (blue) regions 20 residues or longer and labeled P at the percent_cutoff or higher
      i=1
      count_regions=0
20    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'P') count_p=count_p+1
      enddo
30    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'P') then
            count_p=count_p+1
         endif
         goto 30
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         p_start(count_regions)=p_region_start
         p_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 20
      count_regions_p=count_regions

c find ID (red) regions 20 residues or longer and labeled D at the percent_cutoff or higher
      i=1
      count_regions=0
40    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'D') count_p=count_p+1
      enddo
50    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'D') then
            count_p=count_p+1
         endif
         goto 50
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         d_start(count_regions)=p_region_start
         d_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 40
      count_regions_d=count_regions

c find F (black) regions 20 residues or longer and labeled F at the percent_cutoff or higher
      i=1
      count_regions=0
60    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'F') count_p=count_p+1
      enddo
70    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'F') then
            count_p=count_p+1
         endif
         goto 70
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         f_start(count_regions)=p_region_start
         f_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 60
      count_regions_f=count_regions

c find first domain (earliest in sequence) and identify its first residue
      count_domains=count_regions_p+count_regions_d+count_regions_f
      low=1000000
      j=1
      if(j.le.count_domains) then
      do i=1,count_regions_p
         if(p_start(i).lt.low) low=p_start(i)
      enddo
      do i=1,count_regions_d
         if(d_start(i).lt.low) low=d_start(i)
      enddo
      do i=1,count_regions_f
         if(f_start(i).lt.low) low=f_start(i)
      enddo
      domain(j)=low
      endif

c find first residue of each additional domain in successive order
80    j=j+1
      low=1000000
      if(j.le.count_domains) then
      do i=1,count_regions_p
         if(p_start(i).lt.low.and.p_start(i).gt.domain(j-1)) 
     &      low=p_start(i)
      enddo
      do i=1,count_regions_d
         if(d_start(i).lt.low.and.d_start(i).gt.domain(j-1)) 
     &      low=d_start(i)
      enddo
      do i=1,count_regions_f
         if(f_start(i).lt.low.and.f_start(i).gt.domain(j-1)) 
     &      low=f_start(i)
      enddo
      domain(j)=low
      goto 80
      endif

c correct for domain-domain overlap when present
      if(count_domains.gt.1) then
c find overlapping domains and split the overlap, which is percent_cutoff*20/2
      do i=1,count_domains-1
      do j=1,count_regions_p
         if(p_start(j).eq.domain(i)) then
         if(p_end(j).ge.domain(i+1)) then
            p_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_d
            if(d_start(jj).eq.domain(i+1)) then
               d_start(jj)=p_end(j)+1
               domain(i+1)=d_start(jj)
            endif
            enddo
            do jj=1,count_regions_f
            if(f_start(jj).eq.domain(i+1)) then
               f_start(jj)=p_end(j)+1
               domain(i+1)=f_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      do j=1,count_regions_d
         if(d_start(j).eq.domain(i)) then
         if(d_end(j).ge.domain(i+1)) then
            d_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_p
            if(p_start(jj).eq.domain(i+1)) then
               p_start(jj)=d_end(j)+1
               domain(i+1)=p_start(jj)
            endif
            enddo
            do jj=1,count_regions_f
            if(f_start(jj).eq.domain(i+1)) then
               f_start(jj)=d_end(j)+1
               domain(i+1)=f_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      do j=1,count_regions_f
         if(f_start(j).eq.domain(i)) then
         if(f_end(j).ge.domain(i+1)) then
            f_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_d
            if(d_start(jj).eq.domain(i+1)) then
               d_start(jj)=f_end(j)+1
               domain(i+1)=d_start(jj)
            endif
            enddo
            do jj=1,count_regions_p
            if(p_start(jj).eq.domain(i+1)) then
               p_start(jj)=f_end(j)+1
               domain(i+1)=p_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      enddo
      endif

      write(*,*)' '
      write(*,'("Number of PS (blue) regions = ",i6)')count_regions_p
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_p
      write(*,'(i4,7x,i6,7x,i6)') i,p_start(i),p_end(i)
      enddo
      write(*,*)' '
      write(*,'("Number of ID (red) regions = ",i7)')count_regions_d
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_d
      write(*,'(i4,7x,i6,7x,i6)') i,d_start(i),d_end(i)
      enddo
      write(*,*)' '
      write(*,'("Number of F (black) regions = ",i5)')count_regions_f
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_f
      write(*,'(i4,7x,i6,7x,i6)') i,f_start(i),f_end(i)
      enddo

      write(*,*)' '
      do i=1,count_domains
      do j=1,count_regions_p
         if(p_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", PS (blue)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,p_start(j),
     &    p_end(j),(p_end(j)-p_start(j)+1)
         endif
      enddo
      do j=1,count_regions_d
         if(d_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", ID  (red)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,d_start(j),
     &    d_end(j),(d_end(j)-d_start(j)+1)
         endif
      enddo
      do j=1,count_regions_f
         if(f_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", F (black)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,f_start(j),
     &    f_end(j),(f_end(j)-f_start(j)+1)
         endif
      enddo
      enddo

      end
