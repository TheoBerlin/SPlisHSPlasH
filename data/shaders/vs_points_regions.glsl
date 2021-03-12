#version 330

uniform mat4 modelview_matrix;
uniform mat4 projection_matrix;

uniform float radius;
uniform float viewport_width;

layout (location = 0) in vec3 position;
layout (location = 1) in uint level;
layout (location = 2) in uint isBorder;

out block
{
	flat vec4 mv_pos;
	flat uint level;
	flat uint isBorder;
}
Out;

void main()
{
    vec4 mv_pos = modelview_matrix * vec4(position, 1.0);
    vec4 proj = projection_matrix * vec4(radius, 0.0, mv_pos.z, mv_pos.w);
    gl_PointSize = viewport_width * proj.x / proj.w;

	Out.mv_pos = mv_pos;
    Out.level = level;
    Out.isBorder = isBorder;
    gl_Position = projection_matrix * mv_pos;
}
