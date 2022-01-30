#include <algorithm>

#include <optional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <list>

#include <thread>
#include <locale>

#include <GL/freeglut.h>
#include "multidimentional_point.h"

constexpr int list_size = 1024*4;
constexpr int base_window_size_x = 720;
constexpr int base_window_size_y = 720;

bool is_first_boot = true;
bool is_animation_active = true;

int window_size_x = base_window_size_x, window_size_y = base_window_size_y;
float inner_space_range = 200;

inline float random_float(float range)
{
	return ((0.f - range) + ((float)rand() / (RAND_MAX / (2 * range))));
}

inline float random_positive_float(float max)
{
	return (((float)rand() / RAND_MAX) * max);
}

using fptype = double;

struct pendulum {
	using pt4 = dixelu::point<4, fptype>;

	fptype g;
	fptype length, mass;
	fptype theta1, theta2, p1, p2;
	fptype x, y;
	pt4 color;

	int observed_index = 0;

	pt4 get_drivative_at(pt4 point) const {
		fptype o_theta1 = point[0];
		fptype o_theta2 = point[1];
		fptype ptheta1 = point[2];
		fptype ptheta2 = point[3];

		constexpr fptype six = 6;
		constexpr fptype two = 2;
		constexpr fptype eight = 8;
		constexpr auto three = (six / two);
		constexpr auto nine = three * three;
		constexpr auto one = two / two;

		fptype theta1dot = (six / (mass * length * length)) * (two * ptheta1 - three * ptheta2 * std::cos(o_theta1 - o_theta2)) / ((eight * two) - nine * std::pow(std::cos(o_theta1 - o_theta2), two));
		fptype theta2dot = (six / (mass * length * length)) * (eight * ptheta2 - three * ptheta1 * std::cos(o_theta1 - o_theta2)) / ((eight * two) - nine * std::pow(std::cos(o_theta1 - o_theta2), two));

		fptype p1dot = - (one / two) * (mass * length * length) * (theta1dot * theta2dot * std::sin(o_theta1 - o_theta2) + three * (g / length) * std::sin(o_theta1));
		fptype p2dot = - (one / two) * (mass * length * length) * (- one * theta1dot * theta2dot * std::sin(o_theta1 - o_theta2) + (g / length) * std::sin(o_theta2));

		return pt4({ theta1dot, theta2dot, p1dot, p2dot });
	}
	inline pt4 get_params() const {
		return pt4{ theta1, theta2, p1, p2 };
	}
	inline void set_params(pt4 point) {
		theta1 = point[0];
		theta2 = point[1];
		p1 = point[2];
		p2 = point[3];
	}

	inline void evaluate(fptype h = 0.01) {
		auto params = get_params();

		auto k1 = get_drivative_at(params);
		auto k2 = get_drivative_at(params + 0.5 * h * k1);
		auto k3 = get_drivative_at(params + 0.5 * h * k2);
		auto k4 = get_drivative_at(params + h * k3);

		params += (h / 6) * (k1 + 2. * k2 + 2. * k3 + k4);
		set_params(params);
	}

	void draw() {
		double x1 = x + length * std::sin(theta1);
		double x2 = x1 + length * std::sin(theta2);
		double y1 = y - length * std::cos(theta1);
		double y2 = y1 - length * std::cos(theta2);

		glColor4d(color[0], color[1], color[2], color[3]);
		glBegin(GL_LINE_STRIP);
		glVertex2d(x, y);
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glEnd();
	}
};

auto hsva_to_rgba = [](pendulum::pt4 hsva)->pendulum::pt4 {
	fptype hh, p, q, t, ff;
	int i;

	if (hsva[1] <= 0)        // < is bogus, just shuts up warnings
		return { hsva[2], hsva[2], hsva[2], hsva[4] };

	hh = hsva[0];
	if (hh >= 360) hh = 0;
	hh /= 60;
	i = (long)hh;
	ff = hh - i;
	p = hsva[2] * (1 - hsva[1]);
	q = hsva[2] * (1 - (hsva[1] * ff));
	t = hsva[2] * (1 - (hsva[1] * (1 - ff)));

	switch (i) {
	case 0:
		return { hsva[2], t, p, hsva[4] };
	case 1:
		return { q, hsva[2], p, hsva[4] };
	case 2:
		return { p, hsva[2], t, hsva[4] };
	case 3:
		return { p, q, hsva[2], hsva[4] };
	case 4:
		return { t, p, hsva[2], hsva[4] };
	default:
		break;
	}
	return { hsva[2], p, q, hsva[4] };
};

void check_for_interpolation(
	std::list<pendulum>& pendulums,
	std::list<pendulum>::iterator cur_it
) {
	if (pendulums.size() < 5)
		throw std::runtime_error("Less than 5 pendulums");

	fptype angle_eps = 0.01 + random_float(0.0001f );
	auto last_it = std::prev(pendulums.end());
	auto next_it = std::next(cur_it);
	bool is_last = cur_it == last_it;

	if (is_last)
		return;

	auto absdtheta1 = (std::abs)(next_it->theta1 - cur_it->theta1);
	auto absdtheta2 = (std::abs)(next_it->theta2 - cur_it->theta2);

	if (absdtheta1 < angle_eps && absdtheta2 < angle_eps)
		return;

	auto p1 = (cur_it)->get_params();
	auto p2 = (next_it)->get_params();

	fptype max_diff = std::max(std::abs(p1[0] - p2[0]), std::abs(p1[1] - p2[1]));
	auto n_prob = max_diff / angle_eps + 1;
	int n = int(std::floor(n_prob));
	if (random_positive_float(1) < (n_prob - n))
		n++;

	auto simple_interpolation = [&p1, &p2](double x) -> pendulum::pt4 {
		return p2 * x + p1 * (1. - x);
	};

	for (int i = 1; i < n; ++i) {
		fptype alpha = fptype(i) / n;
		pendulum pd = *cur_it;
		auto interpolated = simple_interpolation(alpha);
		pd.set_params(std::move(interpolated));
		pendulums.insert(next_it, pd);
	}

	if (cur_it->observed_index > (pendulums.size() >> 1))
		while (pendulums.size() > list_size)
			pendulums.pop_front();
	else
		while (pendulums.size() > list_size)
			pendulums.pop_back();
}

void evalute_set_of_pendulums(std::list<pendulum>& pendulums) {
	size_t idx = 0;
	for (auto it = pendulums.begin(); it != pendulums.end(); ++it)
		it->evaluate(0.0125);

	for (auto it = pendulums.begin(); it != pendulums.end(); ++it, ++idx) {
		check_for_interpolation(pendulums, it);
	}
}

////////////////////////////////////////////////
///////PULLED FROM SAFC WITHOUT CHANGES/////////
////////////////////////////////////////////////

#define _WH(MainWindow,Element) ((*(*WH)[MainWindow])[Element])//...uh
#define _WH_t(MainWindow,Element,Type) ((Type)_WH(MainWindow,Element))

void onTimer(int v);
void mDisplay() {
	static std::list<pendulum> pends;

	glClear(GL_COLOR_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	auto color_by_index = [](fptype i) -> decltype(pendulum::color) {
		auto t = fptype(i) / list_size;
		return { 0.5f * t * t, 0.5f * (1 + t), 1.f, 0.25f * t * (1.f - t)};
	};

	if (is_first_boot) {
		is_first_boot = false;

		for (int i = 0; i < list_size; i++) {
			pendulum pd;
			const fptype w = 100;
			pd.x = 0;//(0.5 - double(i)/list_size)*2.*w;
			pd.y = 0;
			pd.color = color_by_index(i);
			pd.g = 10;
			pd.length = 75;
			pd.mass = 10;
			pd.p1 = 0;
			pd.p2 = 0;
			pd.theta1 = 3.1;
			pd.theta2 = 2.9 + (1e-1 / list_size) * (fptype(i) + 1);
			pd.observed_index = i;
			pends.push_back(pd);
		}

		//ANIMATION_IS_ACTIVE = !ANIMATION_IS_ACTIVE;
		onTimer(0);
	}

	if(is_animation_active)
		evalute_set_of_pendulums(pends);

	int cnt = 0;
	for (auto& p : pends) {
		p.observed_index = cnt;
		p.color = color_by_index(cnt);
		p.draw();
		cnt++;
	}

	glutSwapBuffers();
}

void mInit() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(
		(0 - inner_space_range) * (window_size_x / base_window_size_x),
		inner_space_range * (window_size_x / base_window_size_x), 
		(0 - inner_space_range) * (window_size_y / base_window_size_y),
		inner_space_range * (window_size_y / base_window_size_y));
}

void onTimer(int v) {
	glutTimerFunc(16, onTimer, 0);
	mDisplay();
}

void OnResize(int x, int y) {
	window_size_x = x;
	window_size_y = y;
	mInit();
	glViewport(0, 0, x, y);
}

void mKey(BYTE k, int x, int y) {
	if (k == 's') 
		is_animation_active = !is_animation_active;

	if (k == 27)
		exit(1);
}
void mSpecialKey(int Key, int x, int y) {
	auto modif = glutGetModifiers();

	if (modif == GLUT_ACTIVE_ALT && Key == GLUT_KEY_DOWN) {
		inner_space_range *= 1.1f;
		OnResize(window_size_x, window_size_y);
	}
	else if (modif == GLUT_ACTIVE_ALT && Key == GLUT_KEY_UP) {
		inner_space_range /= 1.1f;
		OnResize(window_size_x, window_size_y);
	}
}

int main(int argc, char** argv) {
	ShowWindow(GetConsoleWindow(), SW_HIDE);

	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);

	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_MULTISAMPLE);
	glutInitWindowSize(base_window_size_x, base_window_size_y);
	glutCreateWindow("Double pendulum test :)");

	glBlendFunc(GL_SRC_ALPHA, GL_ONE);//_MINUS_SRC_ALPHA
	glEnable(GL_BLEND);

	glEnable(GL_LINE_SMOOTH);//GL_POLYGON_SMOOTH
	glEnable(GL_POINT_SMOOTH);

	glShadeModel(GL_SMOOTH);

	glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);//GL_FASTEST//GL_NICEST
	glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

	glutReshapeFunc(OnResize);
	glutSpecialFunc(mSpecialKey);
	glutKeyboardFunc(mKey);
	glutDisplayFunc(mDisplay);
	mInit();
	glutMainLoop();
	return 0;
}