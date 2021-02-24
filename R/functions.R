# ---- Description --------------------------------------------------------------------------------

# This program defines functions for the analysis.

# ---- Date functions -----------------------------------------------------------------------------

yyyymm_dashes = function(date) {
  as.numeric(substr(date, 1, 4)) * 100 + as.numeric(substr(date, 6, 7))
}

mdist = function(beg_date, end_date) {
  a = as.numeric(beg_date)
  b = as.numeric(end_date)
  y = floor(b / 100) - floor(a / 100)
  m = (b %% 100) - (a %% 100)
  y * 12 + m
}

# ---- Payment functions --------------------------------------------------------------------------

calc_pmt = function(amt, annual_rate, months) {
  ifelse(
    near(annual_rate, 0),
    amt / months,
    amt * ((annual_rate / 12) / (1 - (1 + annual_rate / 12)^(-1 * months)))
  )
}

# ---- Cap and floor values -----------------------------------------------------------------------

cap_vals = function(x, floor, cap) {
  ifelse(
    is.infinite(x),
    NA,
    ifelse(
      x < floor,
      floor,
      ifelse(
        x > cap,
        cap,
        x
      )
    )
  )
}

# ---- State functions ----------------------------------------------------------------------------

us_states = function() {
  c(
    "AL",
    "AK",
    "AZ",
    "AR",
    "CA",
    "CO",
    "CT",
    "DE", 
    "FL",
    "GA",
    "HI",
    "ID",
    "IL",
    "IN",
    "IA",
    "KS",
    "KY",
    "LA",
    "ME",
    "MD",
    "MA",
    "MI",
    "MN",
    "MS",
    "MO",
    "MT",
    "NE",
    "NV",
    "NH",
    "NJ",
    "NM",
    "NY",
    "NC",
    "ND",
    "OH",
    "OK",
    "OR",
    "PA",
    "RI",
    "SC",
    "SD",
    "TN",
    "TX",
    "UT",
    "VT",
    "VA",
    "WA",
    "WV",
    "WI",
    "WY",
    "DC"
  )
}

# ---- Lender functions ---------------------------------------------------------------------------

lender_names = function() {
  tibble::tibble(
    ticker = c(
      "AHFC",
      "Ally Bank",
      "BMW Bank of North America",
      "CBS",
      "CONA",
      "ESB",
      "Fifth Third Bank",
      "Ford Credit",
      "GM FINANCIAL",
      "HCA",
      "MBFS USA LLC",
      "Mechanics Bank",
      "NMAC",
      "SC",
      "TMCC",
      "USAAFSB",
      "VW Credit",
      "WORLD OMNI FINANCIAL CORP"
    ),
    name = c(
      "Honda",
      "Ally Bank",
      "BMW",
      "CarMax",
      "Capital One", 
      "Harley Davidson",
      "Fifth Third",
      "Ford",
      "GM",
      "Hyundai",
      "Mercedes",
      "Mechanics Bank",
      "Nissan",
      "Santander",
      "Toyota",
      "USAA",
      "Volkswagen",
      "World Omni"
    ),
    type = c(
      "Captive",
      "Bank",
      "Captive",
      "Other",
      "Bank",
      "Captive",
      "Bank",
      "Captive",
      "Captive",
      "Captive",
      "Captive",
      "Bank",
      "Captive",
      "Santander",
      "Captive",
      "Bank",
      "Captive",
      "Bank"
    )
  )
}

# ---- Polynomial orders --------------------------------------------------------------------------

form_poly = function(x, d, p) {
  if(p > 1){
    paste0(
      form_poly(x, d, p-1),
      " + ",
      "I(", x, "^", p, ")",
      " + ",
      "I(", d, "*", x, "^", p, ")"
    )
  } else if (p == 1) {
    paste0(
      x,
      " + ",
      "I(", d, "*", x, ")"
    )
  }
}

# ---- Plot functions -----------------------------------------------------------------------------

srd_plot = function(df, y, x, xlab = "", ylab = "", c = 0, smooth = TRUE, p_l = 4, p_r = 4) {
  df[, "y"] = df[, y]
  df[, "x"] = df[, x]
  gg1 = ggplot2::ggplot(
    data = df
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = x,
        y = y
      ),
      color = "darkgray"
    ) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1, 
        size = 12
      ),
      axis.text.y = ggplot2::element_text(
        size = 12
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = NA, color = "black"),
      strip.text.x = ggplot2::element_text(
        size = 14
      )
    ) + 
    ggplot2::labs(
      x = xlab,
      y = ylab
    ) + 
    ggplot2::geom_vline(
      xintercept = c,
      color = "black",
      linetype = "dashed"
    )
  if(smooth == TRUE) {
    gg1 = gg1 + 
      ggplot2::geom_smooth(
        data = dplyr::filter(df, x < 0),
        ggplot2::aes(
          x = x,
          y = y
        ),
        method = "lm", 
        formula = y ~ poly(x, p_l, raw = TRUE),
        color = "black",
        se = FALSE
      ) +
      ggplot2::geom_smooth(
        data = dplyr::filter(df, x >= 0),
        ggplot2::aes(
          x = x,
          y = y
        ),
        method = "lm", 
        formula = y ~ poly(x, p_r, raw = TRUE),
        color = "black",
        se = FALSE
      )
  }
  gg1
}

# ---- Sharp regression discontinuity estimation --------------------------------------------------

form_srd = function(y, x, d, p, fe, cluster) {
  paste0(
    y,
    " ~ ",
    d,
    " + ",
    form_poly(x = x, d = d, p = p),
    " | ", 
    ifelse(
      missing(fe), 
      "0",
      paste0(fe, sep = "", collapse = " + ")
    ), 
    " | 0 | ", 
    ifelse(
      missing(cluster), 
      "0",
      paste0(cluster, sep = "", collapse = " + ")
    )
  )
}
  
srd = function(df, y, x, d, p, h, fe, cluster) {
  df = as.data.frame(df)
  if(missing(h)) {
    df = df[!is.na(df[, y]), ]
  } else {
    df = df[(!is.na(df[, y])) & (abs(df[, x]) <= h), ]
  }
  lfe::felm(
    formula(form_srd(
      y = y,
      x = x,
      d = d,
      p = p,
      fe = fe,
      cluster = cluster
    )),
    data = df
  )
}

srd_matrix = function(df, y, x, d, p_vec, h_vec, fe, cluster) {
  n = 1
  for(h in h_vec) {
    for(p in p_vec) {
      temp_df = srd(
        df = df,
        y = y, 
        x = x,
        d = d,
        p = p,
        h = h,
        fe = fe,
        cluster = cluster
      ) %>%
        broom::tidy(
          robust = TRUE
        ) %>% 
        dplyr::filter(
          term == d
        ) %>%
        dplyr::mutate(
          bandwidth = h,
          polynomial = p
        ) %>% 
        dplyr::select(
          bandwidth,
          polynomial,
          estimate,
          std.error,
          statistic
        ) %>%
        dplyr::rename(
          coef = estimate,
          se = std.error,
          t_stat = statistic
        ) 
      if(n == 1) {
        srd_df = temp_df
      } else {
        srd_df = dplyr::bind_rows(
          srd_df,
          temp_df
        )
      }
      n = n + 1
    }
  }
  srd_df
}
         
form_stacked_srd = function(y, xs, ds, p, fe, cluster) {
  paste0(
    y,
    " ~ ",
    paste0(ds, sep = "", collapse = " + "),
    " + ",
    paste0(form_poly(x = xs, d = ds, p = p), sep = "", collapse = " + "),
    " | ", 
    ifelse(
      missing(fe), 
      "0",
      paste0(fe, sep = "", collapse = " + ")
    ), 
    " | 0 | ", 
    ifelse(
      missing(cluster), 
      "0",
      paste0(cluster, sep = "", collapse = " + ")
    )
  )
} 

stacked_srd = function(df, y, x, d, group, p, h, fe, cluster) {
  if(missing(h)) {
    df = df[!is.na(df[, y]), ]
  } else {
    df = df[(!is.na(df[, y])) & (abs(df[, x]) <= h), ]
  }
  group_vals = unique(as.data.frame(df)[, group])
  for(g in group_vals) {
    df[, paste0(d, "_", gsub(" ", "", g))] = df[, d] * as.numeric(df[, group] == g)
    df[, paste0(x, "_", gsub(" ", "", g))] = df[, x] * as.numeric(df[, group] == g)
  }
  ds = paste0(d, "_", gsub(" ", "", group_vals))
  xs = paste0(x, "_", gsub(" ", "", group_vals))
  lfe::felm(
    formula(form_stacked_srd(
      y = y,
      xs = xs,
      ds = ds,
      p = p,
      fe = fe,
      cluster = cluster
    )),
    data = df
  )
}
        
stacked_srd_matrix = function(df, y, x, d, group, p_vec, h_vec, fe, cluster) {
  group_vals = unique(as.data.frame(df)[, group])
  extracts = paste0(d, "_", gsub(" ", "", group_vals))
  n = 1
  for(h in h_vec) {
    for(p in p_vec) {
      temp_df = stacked_srd(
        df = df,
        y = y, 
        x = x,
        d = d,
        group = group,
        p = p,
        h = h,
        fe = fe,
        cluster = cluster
      ) %>%
        broom::tidy(
          robust = TRUE
        ) %>% 
        dplyr::filter(
          term %in% extracts
        ) %>%
        dplyr::mutate(
          bandwidth = h,
          polynomial = p,
          group = gsub(paste0(d, "_"), "", term)
        ) %>% 
        dplyr::select(
          bandwidth,
          polynomial,
          group,
          estimate,
          std.error,
          statistic
        ) %>%
        dplyr::rename(
          coef = estimate,
          se = std.error,
          t_stat = statistic
        ) 
      if(n == 1) {
        srd_df = temp_df
      } else {
        srd_df = dplyr::bind_rows(
          srd_df,
          temp_df
        )
      }
      n = n + 1
    }
  }
  srd_df
}

np_srd = function(df, y, x) {
  rdrobust::rdrobust(
    y = as.data.frame(df)[, y],
    x = as.data.frame(df)[, x]
  )
}

np_srd_output = function(np) {
  tibble::tibble(
    coef = np$Estimate[, "tau.us"],
    se = np$se["Conventional", ],
    t_stat = np$z["Conventional", ],
    bw = np$bws["h", "left"]
  )
}

np_srd_matrix = function(df, y, x, disc_id) {
  n = 1
  df = as.data.frame(df)
  for(i in sort(unique(df[, disc_id]))) {
    temp_df = np_srd_output(
      np_srd(
        df = df[df[, disc_id] == i, ],
        y = y,
        x = x
      )
    )
    if(n == 1) {
      np_srd_df = temp_df
    } else {
      np_srd_df = dplyr::bind_rows(
        np_srd_df,
        temp_df
      ) 
    }
    n = n + 1
  }
  dplyr::bind_cols(
    np_srd_df,
    tibble::tibble(
      disc_id = sort(unique(df[, disc_id]))
    )
  )
}

# ---- Fuzzy regression discontinuity estimation --------------------------------------------------

form_frd = function(y, v, x, d, p, fe, cluster, ctrls) {
  paste0(
    y,
    " ~ ",
    form_poly(x = x, d = d, p = p),
    ifelse(
      missing(ctrls),
      "",
      paste0(" + ", paste(ctrls, sep = "", collapse = " + "))
    ),
    " | ", 
    ifelse(
      missing(fe), 
      "0",
      paste0(fe, sep = "", collapse = " + ")
    ), 
    " | ",
    "(", v, " ~ ", d, ")",
    " | ", 
    ifelse(
      missing(cluster), 
      "0",
      paste0(cluster, sep = "", collapse = " + ")
    )
  )
}

frd = function(df, y, v, x, d, p, h, fe, cluster, ctrls) {
  df = as.data.frame(df)
  if(missing(h)) {
    df = df[!is.na(df[, y]), ]
  } else {
    df = df[(!is.na(df[, y])) & (abs(df[, x]) <= h), ]
  }
  lfe::felm(
    formula(form_frd(
      y = y,
      v = v,
      x = x,
      d = d,
      p = p,
      fe = fe,
      cluster = cluster,
      ctrls = ctrls
    )),
    data = df
  )
}

frd_matrix = function(df, y, v, x, d, p_vec, h_vec, fe, cluster, ctrls) {
  n = 1
  for(h in h_vec) {
    for(p in p_vec) {
      temp_df = frd(
        df = df,
        y = y, 
        v = v,
        x = x,
        d = d,
        p = p,
        h = h,
        fe = fe,
        cluster = cluster,
        ctrls = ctrls
      ) %>%
        broom::tidy(
          robust = TRUE
        ) %>% 
        dplyr::filter(
          term == paste0("`", v, "(fit)`")
        ) %>%
        dplyr::mutate(
          bandwidth = h,
          polynomial = p
        ) %>% 
        dplyr::select(
          bandwidth,
          polynomial,
          estimate,
          std.error,
          statistic
        ) %>%
        dplyr::rename(
          coef = estimate,
          se = std.error,
          t_stat = statistic
        ) 
      if(n == 1) {
        srd_df = temp_df
      } else {
        srd_df = dplyr::bind_rows(
          srd_df,
          temp_df
        )
      }
      n = n + 1
    }
  }
  srd_df
}


form_frd_cs = function(y, v, x, d, cs, p, fe, cluster) {
  paste0(
    y,
    " ~ ",
    form_poly(x = x, d = d, p = p),
    " + ",
    cs,
    " | ", 
    ifelse(
      missing(fe), 
      "0",
      paste0(fe, sep = "", collapse = " + ")
    ), 
    " | ",
    "(", v, " | I(", cs, " * ", v, ")", " ~ ", d, " + I(", d, "*", cs, "))",
    " | ", 
    ifelse(
      missing(cluster), 
      "0",
      paste0(cluster, sep = "", collapse = " + ")
    )
  )
}

frd_cs = function(df, y, v, x, d, cs, p, h, fe, cluster) {
  df = as.data.frame(df)
  if(missing(h)) {
    df = df[!is.na(df[, y]), ]
  } else {
    df = df[(!is.na(df[, y])) & (abs(df[, x]) <= h), ]
  }
  lfe::felm(
    formula(form_frd_cs(
      y = y,
      v = v,
      x = x,
      d = d,
      cs = cs,
      p = p,
      fe = fe,
      cluster = cluster
    )),
    data = df
  )
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    