---
title: "A Differences-in-Differences Approach using Integer Programming \U0001f615"
description: "Diff-in-Diff using Matching"
lead: "Difference-in-differences (DD) is a commonly used approach in policy evaluation for identifying the impact of an intervention or treatment. Under a parallel trend assumption, we can recover a causal effect by comparing the difference in outcomes between a treatment and a control group, both before and after an intervention was set in place. However, time-varying confounders often break the identifying assumption, biasing our estimates. In this paper, I identify the different contexts in which matching can help reduce such biases, and show how balancing covariates directly can yield better results for solving these issues and bound causal estimates. I also show how this method can be applied both for panel and repeated cross-sectional data. I illustrate these results with simulations and a case study of the impact of a new voucher scheme on socioeconomic segregation in Chile."
date: 2020-12-02T09:19:42+01:00
lastmod: 2020-12-02T09:19:42+01:00
draft: false
weight: 50
images: ["say-hello-to-doks.png"]
contributors: ["Magdalena Bennett"]
---

<style>
.resp-container {
    position: relative;
    overflow: hidden;
    padding-top: 56.25%;
}

.testiframe {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    border: 0;
}
</style>

<div class="resp-container">
    <iframe class="testiframe" src="https://maibennettslides.netlify.app/sds_20201016/mbennett_did.html">
      Oops! Your browser doesn't suppor this.
    </iframe>
</div>