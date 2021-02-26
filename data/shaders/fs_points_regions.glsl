#version 330

uniform float radius;
uniform mat4 projection_matrix;

in block
{
	flat vec4 mv_pos;
	flat uint level;
}
In;

out vec4 out_color;

const vec3 levelColors[3] =
vec3[3](
    vec3(1.0, 0.0, 0.0),
    vec3(1.0, 1.0, 0.0),
    vec3(0.0, 0.0, 1.0)
);

void main()
{
    // calculate normal
    vec3 n;
    n.xy = gl_PointCoord* 2.0 - vec2(1.0);
    float mag = dot(n.xy, n.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
    n.z = sqrt(1.0 - mag);

	vec3 eye = In.mv_pos.xyz + vec3(0.0, 0.0, radius * n.z);
	float depth = (projection_matrix[2][2] * eye.z + projection_matrix[3][2])
        / (projection_matrix[2][3] * eye.z + projection_matrix[3][3]);

    gl_FragDepth = (depth + 1.0) / 2.0;

	// compute final color
    out_color = vec4(levelColors[In.level], 1.0);
}
