varying vec4 color;

void main(void)
{
	gl_FragColor = color;
#ifdef ZHACK
	SetFragDepth(gl_TexCoord[0].z);
#endif
}
