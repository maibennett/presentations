---
output:
  xaringan::moon_reader:
    self_contained: true
    css: ["xaringan-themer.css", "style.css"]
    lib_dir: libs
    nature:
      slideNumberFormat: |
        <div class="progress-bar-container">
          <div class="progress-bar" style="width: calc(%current% / %total% * 100%);">
          </div>
        </div>
      highlightStyle: github
      highlightLines: true
      ratio: "16:9"
      beforeInit: "macros.js"
    includes:
      in_header: header.html
    seal: false
always_allow_html: true
tables: true
---
