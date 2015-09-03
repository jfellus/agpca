/*
 * jpg.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef JPG_H_
#define JPG_H_

#include <jpeglib.h>

int COLORMAP_SIZE = 6;
float COLORMAP[][4] {
		{-1, 0,255,0},
		{0,  0,0,0},
		{0.01,  0,0,255},
		{1, 255,0,0},
		{100, 255,255,255}
};


void interpolate_color(float value, unsigned char* rgb) {
	int col=0;
	for(; col<COLORMAP_SIZE-1; col++) {
		if(value<COLORMAP[col+1][0]) break;
	}

	value = (value-COLORMAP[col][0])/(COLORMAP[col+1][0]-COLORMAP[col][0]);
	rgb[0] = (unsigned char) (COLORMAP[col][1]*(1-value) + COLORMAP[col+1][1]*value);
	rgb[1] = (unsigned char) (COLORMAP[col][2]*(1-value) + COLORMAP[col+1][2]*value);
	rgb[2] = (unsigned char) (COLORMAP[col][3]*(1-value) + COLORMAP[col+1][3]*value);
}

void write_jpeg(FILE* f, float* buffer, size_t w, size_t h, float min=0, float max=1) {
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
	int row_stride;		/* physical row width in image buffer */

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	jpeg_stdio_dest(&cinfo, f);

	cinfo.image_width = w; 	/* image width and height, in pixels */
	cinfo.image_height = h;
	cinfo.input_components = 3;		/* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 80, TRUE /* limit to baseline-JPEG values */);

	jpeg_start_compress(&cinfo, TRUE);
	row_stride = w * 3;

	row_pointer[0] = new JSAMPLE[row_stride];
	unsigned char c[3];
	while (cinfo.next_scanline < cinfo.image_height) {
		for(int i=0; i<w; i++) {
			float v = buffer[cinfo.next_scanline * w + i];
			/*float v = (buffer[cinfo.next_scanline * w + i]-min)/(max-min);
			if(v<0) v=0;
			if(v>1) v=1;*/
			interpolate_color(v, c);
			row_pointer[0][i*3] = c[0];
			row_pointer[0][i*3+1] = c[1];
			row_pointer[0][i*3+2] = c[2];
		}
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	delete row_pointer[0];

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
}



#endif /* JPG_H_ */
