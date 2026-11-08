# Nithin Venkatramanan, personal site

Static site with no build step. To preview locally:

```sh
python3 -m http.server 8000   # then open http://localhost:8000
```

## Files
- `index.html`: the home page
- `tissue.html`, `spectra.html`, `exclusivity.html`: the three interactive figure pages
- `write.html`: the writings list
- `photos.html`: the photo carousel, grid, and full-size view; web-sized images are in `assets/photos/`
- `assets/site.css`: colors and type are set at the top in `:root`
- `assets/figures.js`: the home page figure, the three interactive figures, the rotating contact link, the bio notes, the writings list, the photo carousel, the back link on inner pages, and the footer clock

## Editing
To add an entry, copy a `<div class="entry">` block in `index.html` and change the text.

To update the résumé, replace the PDF and regenerate `assets/resume.png` from it (the pop-up shows the image; Download gives the PDF).

To add a day of writing, copy an `<li class="day">` block to the top of the list in `write.html`. Numbers like 01 / are added automatically. Change the date, title, text, and the id, which appears twice and must be unique.

## Deploying
Any static host works, for example GitHub Pages, Netlify, Vercel, or Cloudflare Pages.
