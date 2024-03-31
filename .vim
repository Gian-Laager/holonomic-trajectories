nnoremap <leader><leader>R :T cargo run --release<CR>
nnoremap <leader><leader>r :T cargo run<CR>
nnoremap <leader><leader>d :T cargo build 2>> /dev/null<CR>:DapContinue<CR>
nnoremap <leader><leader>t :T cargo test<CR>
nnoremap <leader><leader>T :T cargo test --release<CR>
nnoremap <leader><leader>b :T cargo build<CR>
nnoremap <leader><leader>B :T cargo build --release<CR>

lua <<EOF
require("dap").adapters.lldb = {
	type = "executable",
	command = "/usr/bin/lldb-vscode", -- adjust as needed
	name = "lldb",
}

local lldb = {
	name = "Launch lldb",
	type = "lldb", -- matches the adapter
	request = "launch", -- could also attach to a currently running process
	program = "target/debug/ray_tracing",
	cwd = "${workspaceFolder}",
	stopOnEntry = false,
	args = {},
	runInTerminal = false,
}

require('dap').configurations.rust = {
	lldb -- different debuggers or more configurations can be used here
}
EOF

nnoremap <leader>b :DapToggleBreakpoint<CR>
