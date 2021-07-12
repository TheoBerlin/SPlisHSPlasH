#version 330

uniform float radius;
uniform mat4 projection_matrix;

in block
{
	flat vec4 mv_pos;
	flat uint level;
	flat uint isBorder;
}
In;

out vec4 out_color;

const uint NUM_COLORS = 4u;

const vec3 levelColors[NUM_COLORS] =
vec3[NUM_COLORS](
    vec3(1.0, 0.0, 0.0),
    vec3(1.0, 1.0, 0.0),
    vec3(0.0, 0.0, 1.0),
    vec3(0.3, 0.3, 0.3) // Region border color
);

void main()
{
    // calculate normal
    vec3 n;
    n.xy = gl_PointCoord* 2.0 - vec2(1.0);
    float mag = dot(n.xy, n.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
    n.z = sqrt(1.0 - mag);

    // calculate lighting
	const vec3 light_dir = vec3(0.0, 0.0, 1.0);
    float diffuse = max(0.0, dot(light_dir, n));

    vec3 eye = In.mv_pos.xyz + vec3(0.0, 0.0, radius * n.z);
    vec3 halfVector = normalize( eye + light_dir);
    float spec = pow(max(0.0, dot(n,halfVector)), 100.0);

	float depth = (projection_matrix[2][2] * eye.z + projection_matrix[3][2])
        / (projection_matrix[2][3] * eye.z + projection_matrix[3][3]);

    gl_FragDepth = (depth + 1.0) / 2.0;

	// compute final color
    uint colorIndex = In.level * (1u - In.isBorder) + (NUM_COLORS - 1u) * In.isBorder;
    vec3 fluidColor = levelColors[colorIndex];

    // compute final color
	vec3 color_ = 0.25 * fluidColor;
    color_ += 0.7 * diffuse * fluidColor;
    color_ += 0.05 * spec * vec3(1.0);
	color_ = clamp(color_, 0.0, 1.0);

    out_color = vec4(color_, 1.0);
}
