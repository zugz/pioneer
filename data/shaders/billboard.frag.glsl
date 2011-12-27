uniform sampler2D some_texture;

void main(void)
{
	gl_FragColor = gl_Color * texture2D(some_texture, gl_TexCoord[1].st);

#ifdef ZHACK
	SetFragDepth(gl_TexCoord[0].z);
#endif
}
