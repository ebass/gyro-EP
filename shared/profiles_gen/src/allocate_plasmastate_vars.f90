subroutine allocate_plasmastate_vars

  use prgen_globals

  implicit none

  allocate(plst_alla_name(plst_dp1_nspec_alla))
  allocate(plst_all_name(plst_dp1_nspec_all))
  allocate(plst_q_all(plst_dp1_nspec_all))
  allocate(plst_m_all(plst_dp1_nspec_all))
  allocate(plst_ns(nx,plst_dp1_nspec_all))
  allocate(plst_ts(nx,plst_dp1_nspec_all))
  allocate(plst_nb(nx))
  allocate(plst_nmini(nx))
  allocate(plst_nfusi(nx))
  allocate(plst_tb(nx))
  allocate(plst_tmini(nx))
  allocate(plst_tfusi(nx))
  allocate(plst_epar(nx))
  allocate(plst_eperp(nx))
  allocate(plst_vol(nx))
  allocate(plst_rho(nx))
  allocate(plst_grho1(nx))
  allocate(plst_phit(nx))
  allocate(plst_psipol(nx))
  allocate(plst_elong(nx))
  allocate(plst_triang(nx))
  allocate(plst_iota(nx))
  allocate(plst_r_midp_in(nx))
  allocate(plst_r_midp_out(nx))
  allocate(plst_z_midp(nx))
  allocate(plst_zeff(nx))
  allocate(plst_epot(nx))
  allocate(plst_omegat(nx))

  allocate(plst_pe_trans(nx))
  allocate(plst_pi_trans(nx))
  allocate(plst_peech(nx))
  allocate(plst_pmini(nx))
  allocate(plst_pminth(nx))
  allocate(plst_picth(nx))
  allocate(plst_pmine(nx))
  allocate(plst_pohme(nx))
  allocate(plst_qie(nx))
  allocate(plst_pbe(nx))
  allocate(plst_pbi(nx))
  allocate(plst_pbth(nx))
  allocate(plst_pfusi(nx))
  allocate(plst_pfusth(nx))
  allocate(plst_pfuse(nx))
  allocate(plst_prad_br(nx))
  allocate(plst_prad_cy(nx))
  allocate(plst_prad_li(nx))
  allocate(plst_tq_trans(nx))
  allocate(plst_sn_trans(nx))

end subroutine allocate_plasmastate_vars
